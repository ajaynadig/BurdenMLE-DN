#' Fit a BurdenMLE-DN model
#'
#' Fits a discrete-mixture model to gene-level de novo variant counts using
#' MixSQP. The input contains one row per gene and, at minimum, observed case
#' counts and gene-specific mutation rates.
#'
#' @param input_data A data frame with one row per gene and columns
#'   `case_count`, `case_rate`, and `N`. See
#'   the repository's "Worked example" page.
#' @param features Optional one-hot numeric matrix assigning genes to
#'   annotation strata. Row names and order must match `input_data`.
#' @param component_endpoints Optional numeric vector of upper endpoints for
#'   the uniform log-rate-ratio mixture components. Use
#'   [effect_size_grid()] to construct prevalence-derived or conservatively
#'   wide grids.
#' @param no_cpts Number of mixture components when `component_endpoints` is
#'   not supplied.
#' @param grid_size Number of numerical-integration points per component.
#' @param mutvar_est Compute mutational variance and annotation summaries.
#' @param max_iter,max_iter_boot,tol Controls retained for the optional legacy
#'   optimizer.
#' @param prevalence Population prevalence on the 0--1 scale.
#' @param bootstrap Run a gene-level nonparametric bootstrap.
#' @param bootstrap_samples Optional matrix of user-supplied bootstrap row
#'   indices.
#' @param n_boot Number of bootstrap replicates.
#' @param null_sim Run optional parametric null simulations. These are not
#'   needed for estimation or bootstrap confidence intervals.
#' @param n_null Number of null simulations.
#' @param return_likelihood Store the maximized log likelihood.
#' @param estimate_posteriors Compute gene-level posterior summaries.
#' @param estimate_effective_penetrance Compute effective penetrance.
#' @param optimizer Optimization routine. MixSQP is the supported default.
#' @param mixsqp_control Named list passed to `mixsqp::mixsqp()`. Package
#'   defaults preserve small components and scale the active-set allowance to
#'   the requested component grid.
#'
#' @return An object of class `BurdenMLEDN_fit`. Important fields include
#'   `delta`, `component_endpoints`, `ll`, `mutvar_output`, `penetrance`,
#'   `bootstrap_output`, and `optimizer_elapsed`.
#' @export
#'
#' @examples
#' example_file <- system.file("extdata", "synthetic_example.csv",
#'                            package = "BurdenMLEDN")
#' dat <- read.csv(example_file, row.names = "gene")
#' fit <- BurdenMLE_DN(
#'   dat,
#'   prevalence = 0.02,
#'   bootstrap = FALSE,
#'   null_sim = FALSE,
#'   estimate_effective_penetrance = TRUE
#' )
#' fit
BurdenMLE_DN <- function(input_data,
                          features = NULL,
                          component_endpoints = NULL,
                          no_cpts = 10,
                          grid_size = 10,
                          mutvar_est = TRUE,
                          max_iter =10000,
                          max_iter_boot = 1000,
                          tol = 1e-6,
                          prevalence = NULL,
                          bootstrap = TRUE,
                          bootstrap_samples = NULL,
                          n_boot = 100,
                          null_sim = FALSE,
                          n_null = 100,
                          return_likelihood = TRUE,
                          estimate_posteriors = TRUE,
                          estimate_effective_penetrance = TRUE,
                          optimizer = c("mixsqp", "EM"),
                          mixsqp_control = list()) {

  optimizer <- match.arg(optimizer)

  if (is.null(prevalence) || length(prevalence) != 1L ||
      !is.finite(prevalence) || prevalence <= 0 || prevalence >= 1) {
    stop("prevalence must be one finite number strictly between 0 and 1.")
  }
  if (length(no_cpts) != 1L || no_cpts < 2 ||
      no_cpts != as.integer(no_cpts)) {
    stop("no_cpts must be one integer of at least 2.")
  }
  if (length(grid_size) != 1L || grid_size < 1 ||
      grid_size != as.integer(grid_size)) {
    stop("grid_size must be one positive integer.")
  }

  genetic_data = process_data_trio(input_data,
                                   features)

  component_endpoints = choose_component_endpoints_trio(component_endpoints,
                                                        no_cpts,
                                                        prevalence)
  cat("...initializing model")
  model = initialize_model(likelihood_function = poisson_uniform_likelihood,
                           genetic_data = genetic_data,
                           component_endpoints = component_endpoints,
                           features = features,
                           grid_size = grid_size)
  model$optimizer <- optimizer


  cat("...running ", optimizer, " in full dataset", sep = "")
  fit_start <- proc.time()[["elapsed"]]
  if (optimizer == "EM") {
    model = EM_fit(model,
                   max_iter,
                   tol = tol,
                   return_likelihood = return_likelihood)
  } else {
    model = MixSQP_fit(model,
                       control = mixsqp_control,
                       return_likelihood = return_likelihood)
  }
  model$optimizer_elapsed <- c(
    full_fit = proc.time()[["elapsed"]] - fit_start
  )

  # Bootstrap the selected optimizer.
  if (bootstrap) {
    bootstrap_start <- proc.time()[["elapsed"]]
    if (optimizer == "EM") {
      model$bootstrap_output = bootstrap_EM(model,
                                            n_boot,
                                            max_iter_boot,
                                            bootstrap_samples,
                                            tol = tol)
    } else {
      model$bootstrap_output = bootstrap_MixSQP(model,
                                                n_boot,
                                                bootstrap_samples,
                                                control = mixsqp_control)
    }
    model$optimizer_elapsed["bootstrap"] <-
      proc.time()[["elapsed"]] - bootstrap_start

  }

  if (null_sim) {
    if (optimizer == "EM") {
      model$null_delta = null_EM_trio(genetic_data,
                                      model,
                                      max_iter,
                                      n_null,
                                      grid_size,
                                      tol)
    } else {
      model$null_delta = null_MixSQP_trio(genetic_data,
                                          model,
                                          n_null,
                                          grid_size,
                                          control = mixsqp_control)
    }
  }

  if (mutvar_est) {
    #estimate mutvar in the full dataset
    model$mutvar_output = estimate_mutvar_trio(model = model,
                                                           genetic_data = genetic_data,
                                                           prevalence = prevalence)

    #bootstrap mutvar estimation
    if (bootstrap) {
      bootstrap_mutvar_output <- lapply(1:n_boot,
                                              function(iter) {



                                                model_boot = model
                                                model_boot$conditional_likelihood = model_boot$conditional_likelihood[model$bootstrap_output$bootstrap_samples[,iter],]
                                                model_boot$features = model_boot$features[model$bootstrap_output$bootstrap_samples[,iter],]
                                                model_boot$delta = model$bootstrap_output$bootstrap_delta[[iter]]

                                                boot_mutvar = estimate_mutvar_trio(model = model_boot,
                                                                                               genetic_data = genetic_data[model$bootstrap_output$bootstrap_samples[,iter],],
                                                                                               prevalence = prevalence)

                                              })

      bootstrap_mutvar_ests = sapply(1:length(bootstrap_mutvar_output), function(x) bootstrap_mutvar_output[[x]]$total_mutvar)
      mutvar_CI = quantile(bootstrap_mutvar_ests,c(0.025,0.975))

      bootstrap_annotmutvar_ests = sapply(1:length(bootstrap_mutvar_output), function(x) bootstrap_mutvar_output[[x]]$annot_mutvar)
      annot_mutvar_CI = sapply(1:nrow(bootstrap_annotmutvar_ests),
                           function(i) {
                             quantile(bootstrap_annotmutvar_ests[i,],c(0.025,0.975))
                           })


      bootstrap_fracmutvar_ests = sapply(1:length(bootstrap_mutvar_output), function(x) bootstrap_mutvar_output[[x]]$frac_mutvar)
      bootstrap_enrich_ests = sapply(1:length(bootstrap_mutvar_output), function(x) bootstrap_mutvar_output[[x]]$enrichment)
      normalized_bootstrap_defined <- all(is.finite(bootstrap_fracmutvar_ests)) &&
        all(is.finite(bootstrap_enrich_ests))
      if (normalized_bootstrap_defined) {
        fracmutvar_CI = sapply(1:nrow(bootstrap_fracmutvar_ests),
                              function(i) {
                                quantile(bootstrap_fracmutvar_ests[i,],c(0.025,0.975))
                              })
        enrich_CI = sapply(1:nrow(bootstrap_enrich_ests),
                           function(i) {
                             quantile(bootstrap_enrich_ests[i,],c(0.025,0.975))
                           })
      } else {
        warning(
          "At least one bootstrap replicate has zero total mutational variance; ",
          "fractional mutational variance and enrichment CIs are undefined and ",
          "will not be calculated."
        )
        fracmutvar_CI = matrix(
          NA_real_, nrow = 2, ncol = nrow(bootstrap_fracmutvar_ests),
          dimnames = list(c("2.5%", "97.5%"), rownames(bootstrap_fracmutvar_ests))
        )
        enrich_CI = matrix(
          NA_real_, nrow = 2, ncol = nrow(bootstrap_enrich_ests),
          dimnames = list(c("2.5%", "97.5%"), rownames(bootstrap_enrich_ests))
        )
      }
      model$mutvar_output$mutvar_CI = mutvar_CI
      model$mutvar_output$annot_mutvar_CI = annot_mutvar_CI
      model$mutvar_output$fracmutvar_CI = fracmutvar_CI
      model$mutvar_output$enrich_CI = enrich_CI

    }

    #null mutvar estimates
    if (null_sim) {
      null_mutvar_ests <- sapply(1:n_null,
                             function(iter) {
                               model_null = model
                               model_null$delta = model$null_delta[[iter]]

                               null_mutvar = estimate_mutvar_trio(model = model_null,
                                                                              genetic_data = genetic_data,
                                                                              prevalence = prevalence)

                               return(null_mutvar$total_mutvar)

                             })

      model$mutvar_output$null_mutvar_ests = null_mutvar_ests
      model$mutvar_output$mutvar_p = mean(
        model$mutvar_output$null_mutvar_ests > model$mutvar_output$total_mutvar
      )
    }


  }

  #Get some posterior expectations and hypothesis tests
  if (estimate_posteriors == TRUE) {
    #first, get the naive rate ratio estimates
    RR_naive = genetic_data$case_count / genetic_data$expected_count

    #get some simple poisson p values
    RR_poisson_p = ppois(genetic_data$case_count,
                         genetic_data$expected_count,
                         lower.tail = FALSE)

    #get posterior means
    RR_posterior_means <- posterior_expectation(model,
                                                genetic_data,
                                                exp,
                                                grid_size)

    #hypothesis test
    RR_problessthanequal0 <- posterior_expectation(model,
                                                   genetic_data,
                                                   function(x) {x <= 0},
                                                   grid_size)

    posterior_gene_estimate_df <- data.frame(Case_Count = genetic_data$case_count,
                                             Expected_Count = genetic_data$expected_count,
                                             Estimate = RR_naive,
                                             Estimate_Poisson_P = RR_poisson_p,
                                             Posterior_Mean = RR_posterior_means,
                                             Posterior_ProbLessEqualZero = RR_problessthanequal0)

    rownames(posterior_gene_estimate_df) <- rownames(genetic_data)

    model$posterior_gene_estimates = posterior_gene_estimate_df

  }

  if (estimate_effective_penetrance) {
    cat("...computing effective penetrance")

    penetrance_moments = precompute_effective_penetrance_moments(
      model,
      genetic_data
    )
    peneff = effective_penetrance_func(model,
                                       genetic_data,
                                       prevalence,
                                       moments = penetrance_moments)
    peneff_CI = NA
    if (bootstrap) {
      cat("...bootstrap effective penetrance")
      bootstrap_peneff_ests = pbsapply(1:length(model$bootstrap_output$bootstrap_delta),
                                     function(iter) {
                                       model_boot = model
                                       model_boot$conditional_likelihood = model_boot$conditional_likelihood[model$bootstrap_output$bootstrap_samples[,iter],]
                                       model_boot$features = model_boot$features[model$bootstrap_output$bootstrap_samples[,iter],]
                                       model_boot$delta = model$bootstrap_output$bootstrap_delta[[iter]]

                                       sample_indices = model$bootstrap_output$bootstrap_samples[,iter]
                                       moments_boot = list(
                                         numerator = penetrance_moments$numerator[sample_indices, , drop = FALSE],
                                         denominator = penetrance_moments$denominator[sample_indices, , drop = FALSE]
                                       )

                                       effective_penetrance_func(model_boot,
                                                                 genetic_data[sample_indices,],
                                                                 prevalence,
                                                                 moments = moments_boot)
                                     })

      if (all(is.finite(bootstrap_peneff_ests))) {
        peneff_CI = quantile(bootstrap_peneff_ests,c(0.025,0.975))
      } else {
        warning(
          "At least one bootstrap replicate has undefined effective penetrance; ",
          "the effective-penetrance CI will not be calculated."
        )
        peneff_CI = c("2.5%" = NA_real_, "97.5%" = NA_real_)
      }
    }

    model$penetrance = list(effective_penetrance = peneff,
                            effective_penetrance_CI = peneff_CI)
  }

  model$call <- match.call()
  model$prevalence <- prevalence
  model$input_summary <- list(
    genes = nrow(genetic_data),
    observed_variants = sum(genetic_data$case_count),
    sample_size = unique(genetic_data$N)
  )
  class(model) <- c("BurdenMLEDN_fit", "list")
  model

}
