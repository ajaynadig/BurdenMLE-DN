inclusive_finite_null_pvalue <- function(null_values, observed_value) {
  (1 + sum(null_values >= observed_value)) / (length(null_values) + 1)
}

#' Fit a BurdenMLE-DN model
#'
#' Fits a discrete-mixture model to gene-level de novo variant counts using
#' MixSQP. The input contains one row per gene, observed case counts, and either
#' exact expected counts or the mutation-rate inputs needed to derive them.
#'
#' @param input_data A data frame with at least two genes; explicit, unique,
#'   nonempty, non-missing gene row names; and either columns `case_count` and
#'   `expected_count`, or columns `case_count`, `case_rate`, and `N`. Supplied
#'   expected counts are retained exactly. See the repository's "Worked
#'   example" page.
#' @param features Optional finite one-hot numeric matrix assigning genes to
#'   nonempty annotation strata. It must have explicit, unique, nonempty,
#'   non-missing gene row names identical to `input_data` in the same order and
#'   unique, nonempty, non-missing stratum column names. Reordered feature rows
#'   fail rather than being realigned.
#' @param component_endpoints Optional numeric vector of upper endpoints for
#'   the uniform log-rate-ratio mixture components. Use
#'   [effect_size_grid()] to construct prevalence-derived or conservatively
#'   wide grids. When supplied, these endpoints take precedence and `no_cpts`
#'   is ignored.
#' @param no_cpts Number of mixture components when `component_endpoints` is
#'   not supplied.
#' @param grid_size Number of numerical-integration points per component.
#' @param mutvar_est Compute mutational variance and annotation summaries.
#'   This uses deterministic analytic moments of the continuous-uniform
#'   components and requires `case_rate`; expected-count-only inputs must set
#'   this to `FALSE`.
#' @param max_iter,max_iter_boot,tol Controls retained for the optional legacy
#'   optimizer.
#' @param prevalence Population prevalence on the 0--1 scale. Required when
#'   component endpoints are generated, mutational variance is requested, or
#'   effective penetrance is requested; otherwise it may be `NULL`.
#' @param bootstrap Run a gene-level nonparametric bootstrap.
#' @param bootstrap_samples Optional character matrix of sampled gene IDs, with
#'   one input-gene draw per row and one replicate per column. Repeated IDs are
#'   allowed. Numeric row positions are not accepted.
#' @param n_boot Number of bootstrap replicates.
#' @param bootstrap_seed Positive integer seed for generated bootstrap samples.
#'   The caller's RNG state is restored. Ignored when `bootstrap_samples` is
#'   supplied.
#' @param null_sim Run optional parametric null simulations. These are not
#'   needed for estimation or bootstrap confidence intervals.
#' @param n_null Number of null simulations.
#' @param null_seed Positive integer seed for null simulations. The caller's
#'   RNG state is restored.
#' @param return_likelihood Store the maximized log likelihood.
#' @param estimate_posteriors Compute gene-level posterior summaries.
#' @param estimate_effective_penetrance Compute gene-average effective
#'   penetrance. The distinct mutation-weighted estimand is available through
#'   [mutation_weighted_effective_penetrance()].
#' @param optimizer Optimization routine. MixSQP is the supported default.
#' @param mixsqp_control Named list passed to `mixsqp::mixsqp()`. Package
#'   defaults preserve small components and scale the active-set allowance to
#'   the requested component grid.
#'
#' @details Gene-level Poisson p-values are inclusive one-sided probabilities,
#'   `P(X >= observed)`. When null simulation and mutational variance are both
#'   requested, the finite-null p-value includes ties and uses a plus-one
#'   correction, so its minimum is `1 / (n_null + 1)`.
#'
#' @return An object of class `BurdenMLEDN_fit`. Important fields include
#'   `delta`, `component_endpoints`, `ll`, `fit_status`,
#'   `uncertainty_reliable`, `mutvar_output`, `penetrance`,
#'   `bootstrap_output`, `null_output`, `input_summary`, and
#'   `optimizer_elapsed`.
#'   `input_summary$sample_size` is `NA` when `N` was not supplied, and
#'   `input_summary$case_rate_available` reports whether mutation-rate-dependent
#'   estimands are supported. Stratum-indexed weights, mutational-variance
#'   summaries, and confidence intervals retain the feature-column names.
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
                          bootstrap_seed = 1L,
                          null_sim = FALSE,
                          n_null = 100,
                          null_seed = 1L,
                          return_likelihood = TRUE,
                          estimate_posteriors = TRUE,
                          estimate_effective_penetrance = TRUE,
                          optimizer = c("mixsqp", "EM"),
                          mixsqp_control = list()) {

  optimizer <- match.arg(optimizer)

  if (bootstrap || null_sim) {
    had_rng_state <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_rng_state) caller_rng_state <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_rng_state) {
        assign(".Random.seed", caller_rng_state, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
  }

  genetic_data = process_data_trio(input_data)
  if (mutvar_est) {
    require_case_rate_for_mutvar(genetic_data)
  }

  prevalence_required <- is.null(component_endpoints) || mutvar_est ||
    estimate_effective_penetrance
  if (prevalence_required || !is.null(prevalence)) {
    if (is.null(prevalence) || length(prevalence) != 1L ||
        !is.numeric(prevalence) || !is.finite(prevalence) ||
        prevalence <= 0 || prevalence >= 1) {
      stop(
        "prevalence must be one finite number strictly between 0 and 1 when ",
        "component endpoints are generated, mutational variance is requested, ",
        "or effective penetrance is requested."
      )
    }
  }
  if (length(grid_size) != 1L || grid_size < 1 ||
      grid_size != as.integer(grid_size)) {
    stop("grid_size must be one positive integer.")
  }

  features = validate_features_trio(genetic_data, features)

  bootstrap_contract <- NULL
  if (bootstrap) {
    bootstrap_contract <- prepare_bootstrap_samples(
      rownames(genetic_data), n_boot, bootstrap_samples, bootstrap_seed
    )
    n_boot <- length(bootstrap_contract$replicate_ids)
  }
  null_contract <- NULL
  if (null_sim) {
    n_null <- validate_positive_integer(n_null, "n_null")
    null_seeds <- resampling_seeds(null_seed, n_null, "null_seed", "null_")
    null_contract <- list(
      replicate_ids = names(null_seeds), seeds = null_seeds
    )
  }

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
  uncertainty_checks <- logical(0)
  if (bootstrap) {
    bootstrap_start <- proc.time()[["elapsed"]]
    if (optimizer == "EM") {
      bootstrap_fit_status <- bootstrap_EM(
        model,
        bootstrap_contract$indices,
        bootstrap_contract$replicate_ids,
        max_iter_boot,
        tol = tol
      )
    } else {
      bootstrap_fit_status <- bootstrap_MixSQP(
        model,
        bootstrap_contract$indices,
        bootstrap_contract$replicate_ids,
        control = mixsqp_control
      )
    }
    bootstrap_reliable <- assess_uncertainty_fits(
      bootstrap_fit_status, "Bootstrap"
    )
    model$bootstrap_output <- list(
      bootstrap_delta = lapply(bootstrap_fit_status, `[[`, "weights"),
      bootstrap_samples = bootstrap_contract$samples,
      bootstrap_indices = bootstrap_contract$indices,
      fit_status = bootstrap_fit_status,
      replicate_ids = bootstrap_contract$replicate_ids,
      replicate_seeds = bootstrap_contract$seeds,
      sample_source = bootstrap_contract$source,
      reliable = bootstrap_reliable
    )
    uncertainty_checks <- c(uncertainty_checks, bootstrap_reliable)
    model$optimizer_elapsed["bootstrap"] <-
      proc.time()[["elapsed"]] - bootstrap_start

  }

  if (null_sim) {
    if (optimizer == "EM") {
      null_fit_status <- null_EM_trio(
        genetic_data,
        model,
        max_iter,
        null_contract$seeds,
        null_contract$replicate_ids,
        grid_size,
        tol
      )
    } else {
      null_fit_status <- null_MixSQP_trio(
        genetic_data,
        model,
        null_contract$seeds,
        null_contract$replicate_ids,
        grid_size,
        control = mixsqp_control
      )
    }
    null_reliable <- assess_uncertainty_fits(null_fit_status, "Null simulation")
    model$null_delta <- lapply(null_fit_status, `[[`, "weights")
    model$null_output <- list(
      null_delta = model$null_delta,
      fit_status = null_fit_status,
      replicate_ids = null_contract$replicate_ids,
      replicate_seeds = null_contract$seeds,
      reliable = null_reliable
    )
    uncertainty_checks <- c(uncertainty_checks, null_reliable)
  }
  model$uncertainty_reliable <- if (length(uncertainty_checks) == 0L) {
    NA
  } else {
    all(uncertainty_checks)
  }

  if (mutvar_est) {
    stratum_names <- colnames(model$features)
    #estimate mutvar in the full dataset
    model$mutvar_output = estimate_mutvar_trio(model = model,
                                                           genetic_data = genetic_data,
                                                           prevalence = prevalence)

    #bootstrap mutvar estimation
    if (bootstrap) {
      bootstrap_mutvar_output <- lapply(seq_len(n_boot), function(iter) {
        replicate <- reconstruct_bootstrap_replicate(model, genetic_data, iter)
        estimate_mutvar_trio(
          model = replicate$model,
          genetic_data = replicate$genetic_data,
          prevalence = prevalence
        )
      })
      names(bootstrap_mutvar_output) <- model$bootstrap_output$replicate_ids
      bootstrap_matrix <- function(field) {
        output <- matrix(
          vapply(
            bootstrap_mutvar_output,
            function(value) unname(value[[field]]),
            numeric(length(stratum_names))
          ),
          nrow = length(stratum_names),
          ncol = length(bootstrap_mutvar_output)
        )
        dimnames(output) <- list(
          stratum_names, model$bootstrap_output$replicate_ids
        )
        output
      }
      bootstrap_CI <- function(values) {
        output <- matrix(
          apply(values, 1L, quantile, probs = c(0.025, 0.975)),
          nrow = 2L
        )
        dimnames(output) <- list(c("2.5%", "97.5%"), stratum_names)
        output
      }

      bootstrap_mutvar_ests <- vapply(
        bootstrap_mutvar_output, `[[`, numeric(1), "total_mutvar"
      )
      mutvar_CI = quantile(bootstrap_mutvar_ests,c(0.025,0.975))

      bootstrap_annotmutvar_ests <- bootstrap_matrix("annot_mutvar")
      annot_mutvar_CI <- bootstrap_CI(bootstrap_annotmutvar_ests)

      bootstrap_fracmutvar_ests <- bootstrap_matrix("frac_mutvar")
      bootstrap_enrich_ests <- bootstrap_matrix("enrichment")
      normalized_bootstrap_defined <- all(is.finite(bootstrap_fracmutvar_ests)) &&
        all(is.finite(bootstrap_enrich_ests))
      if (normalized_bootstrap_defined) {
        fracmutvar_CI <- bootstrap_CI(bootstrap_fracmutvar_ests)
        enrich_CI <- bootstrap_CI(bootstrap_enrich_ests)
      } else {
        warning(
          "At least one bootstrap replicate has zero total mutational variance; ",
          "fractional mutational variance and enrichment CIs are undefined and ",
          "will not be calculated."
        )
        fracmutvar_CI = matrix(
          NA_real_, nrow = 2, ncol = nrow(bootstrap_fracmutvar_ests),
          dimnames = list(c("2.5%", "97.5%"), stratum_names)
        )
        enrich_CI = matrix(
          NA_real_, nrow = 2, ncol = nrow(bootstrap_enrich_ests),
          dimnames = list(c("2.5%", "97.5%"), stratum_names)
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
      model$mutvar_output$mutvar_p <- inclusive_finite_null_pvalue(
        model$mutvar_output$null_mutvar_ests,
        model$mutvar_output$total_mutvar
      )
    }


  }

  #Get some posterior expectations and hypothesis tests
  if (estimate_posteriors == TRUE) {
    #first, get the naive rate ratio estimates
    RR_naive = genetic_data$case_count / genetic_data$expected_count

    #get some simple poisson p values
    RR_poisson_p = ppois(genetic_data$case_count - 1,
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
      bootstrap_peneff_ests <- pbsapply(
        seq_along(model$bootstrap_output$bootstrap_delta),
        function(iter) {
          replicate <- reconstruct_bootstrap_replicate(model, genetic_data, iter)
          moments_boot <- list(
            numerator = penetrance_moments$numerator[
              replicate$indices, , drop = FALSE
            ],
            denominator = penetrance_moments$denominator[
              replicate$indices, , drop = FALSE
            ]
          )
          effective_penetrance_func(
            replicate$model,
            replicate$genetic_data,
            prevalence,
            moments = moments_boot
          )
        }
      )
      names(bootstrap_peneff_ests) <- model$bootstrap_output$replicate_ids

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
    sample_size = if ("N" %in% names(genetic_data)) {
      unique(genetic_data$N)
    } else {
      NA_real_
    },
    case_rate_available = "case_rate" %in% names(genetic_data)
  )
  class(model) <- c("BurdenMLEDN_fit", "list")
  model

}
