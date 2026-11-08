draw_dirichlet <- function(alpha) {
  draws <- rgamma(length(alpha), shape = alpha, rate = 1)
  draws / sum(draws)
}

generate_effect_sizes <- function(n_genes, distribution, params) {
  effects <- switch(distribution,
    "null" = rep(0, n_genes),
    "half_uniform" = {
      # Mixture of half-uniforms as used in the model
      component <- sample(1:length(params$endpoints), n_genes, prob = params$weights, replace = TRUE)
      sapply(1:n_genes, function(i) runif(1, min=0, max=params$endpoints[component[i]]))
    },
    "half_infinitesimal" = {
      # Positive half-normal distribution
      abs(rnorm(n_genes, mean=0, sd=params$sd))
    },
    "point_normal" = {
      # Point mass at zero and positive normal
      is_zero <- rbinom(n_genes, 1, params$pi0)
      effects <- abs(rnorm(n_genes, mean=params$mean, sd=params$sd))
      effects[is_zero == 1] <- 0
      effects
    },
    "oligogenic" = {
      # Majority zero, few large effects
      effects <- rep(0, n_genes)
      n_causal <- round(n_genes * params$prop_causal)
      causal_idx <- sample(1:n_genes, n_causal)
      effects[causal_idx] <- params$effect_size  # Single effect size for all causal genes
      effects
    }
  )
  return(effects)
}

# Function to simulate counts based on mutation rates and effect sizes
simulate_counts <- function(mu_rates, effect_sizes, N, empirical_rownames) {
  # Simulate case counts using correct Poisson parameter
  case_counts <- rpois(length(mu_rates), 2 * N * mu_rates * exp(effect_sizes))
  
  # Create the input data frame in the correct format for BurdenMLE-DN.
  sim_input <- data.frame(
    case_count = case_counts,
    case_rate = mu_rates,
    N = rep(N, length(mu_rates))
  )
  

  sim_input$expected_count <- 2 * sim_input$N * sim_input$case_rate
  # Set rownames to match empirical data
  rownames(sim_input) <- empirical_rownames
  
  return(sim_input)
}

# Function to calculate true mutational variance from effect sizes
calculate_true_mutvar <- function(effect_sizes, mu_rates, prevalence) {
  # Calculate mutational variance per gene and sum
  gamma <- exp(effect_sizes)
  mutvar_per_gene <- (prevalence * (gamma - 1)^2 * 2*mu_rates) / (1-prevalence)
  true_mutvar <- sum(mutvar_per_gene)
  return(true_mutvar)
}

# Function to calculate fraction of cases with causal mutations
calculate_true_frac_cases <- function(effect_sizes, sim_data, threshold) {
  # Count mutations in genes with effect size > log(threshold)
  causal_mutations <- sum(sim_data$case_count[effect_sizes > log(threshold)])
  
  # Get N from the simulation data
  N <- sim_data$N[1]  # All N values are the same
  return(causal_mutations/N)
}

# Calculate true effective penetrance based on true effect sizes
calculate_true_effective_penetrance <- function(true_effects, sim_data, prevalence) {
  rate_ratio <- exp(true_effects)
  mutation_weights <- sim_data$case_rate
  peneff_numerator <- sum(
    mutation_weights * (rate_ratio - 1) * rate_ratio
  )
  peneff_denominator <- sum(mutation_weights * (rate_ratio - 1))
  if (!is.finite(peneff_denominator) || peneff_denominator == 0) {
    return(NA_real_)
  }
  prevalence * peneff_numerator / peneff_denominator
}

get_fraccase<- function(model,genetic_data,threshold) {

  fraccases = sum(2*genetic_data$case_rate*posterior_expectation(model,genetic_data,
                                                                                        function(x) {
                                                                                          exp(x) * (exp(x) > threshold)
                                                                                        },
                                                                                        grid_size = 10))
  return(fraccases)
}

# Run simulation for one scenario
run_simulation <- function(distribution, params, empirical_data, empirical_features, N, n_sims=1,
                           optimizer, progress_tick = NULL, job_label = NULL) {
  results <- list()
  n_genes <- nrow(empirical_data)
  mu_rates <- empirical_data$case_rate  # Get base mutation rates
  
  for(i in 1:n_sims) {
    # Sample new parameters for this iteration if they are functions
    iter_params <- params
    if(distribution != "null") {
      for(param_name in names(params)) {
        if(is.function(params[[param_name]])) {
          iter_params[[param_name]] <- params[[param_name]]()
        }
      }
    }
    
    # Generate true effects
    true_effects <- generate_effect_sizes(n_genes, distribution, iter_params)
    
    # Calculate true mutational variance
    true_mutvar <- calculate_true_mutvar(true_effects, mu_rates, iter_params$prevalence)
    
    # Simulate data with specified N
    sim_data <- simulate_counts(mu_rates, true_effects, N, rownames(empirical_data))
    
    # Calculate true fraction of cases with causal mutations
    true_frac_cases <- calculate_true_frac_cases(true_effects, sim_data, iter_params$RR_threshold)

    # Calculate true effective penetrance
    true_peneff <- calculate_true_effective_penetrance(true_effects, sim_data, iter_params$prevalence)

    # Some optimizers use randomized numerical routines internally. Preserve the
    # simulation RNG stream so optimizer choice cannot change later datasets.
    rng_state_after_simulation <- .Random.seed

    # Fit the BurdenMLE-DN model using empirical features.
    fit_timing <- system.time(model <- BurdenMLE_DN(
      input_data = sim_data,
      features = empirical_features,
      component_endpoints = seq(0, log(1/0.01), length.out = 10),  # Back to 10 components
      null_sim = FALSE,
      bootstrap = FALSE,
      prevalence = iter_params$prevalence,
      estimate_posteriors = TRUE,
      estimate_effective_penetrance = TRUE,
      optimizer = optimizer
    ))
    assign(".Random.seed", rng_state_after_simulation, envir = .GlobalEnv)

    final_log_likelihood <- model$fit_status$log_likelihood
    optimizer_status <- model$fit_status$backend_message

    #get the estimated fraction of cases with large effect mutations, using the same threshold
    est_fraccases <- get_fraccase(model, sim_data, iter_params$RR_threshold)

    # Store estimates and compact diagnostics, but not the full fitted model.
    results[[i]] <- list(
      true_effects = true_effects,
      simulation_data_hash = digest::digest(sim_data, algo = "sha256"),
      true_mutvar = true_mutvar,
      true_frac_cases = true_frac_cases,
      true_peneff = true_peneff,
      estimated_mutvar = model$mutvar_output$total_mutvar,
      estimated_frac_cases = est_fraccases,
      estimated_peneff = model$penetrance$effective_penetrance,
      model_weights = model$delta[1,],
      model_delta = model$delta,
      optimizer = optimizer,
      log_likelihood = unname(final_log_likelihood),
      log_likelihood_per_gene = unname(final_log_likelihood) / n_genes,
      total_fit_seconds = unname(fit_timing[["elapsed"]]),
      optimizer_full_fit_seconds = unname(model$optimizer_elapsed[["full_fit"]]),
      optimizer_iterations = model$fit_status$iterations,
      optimizer_status = optimizer_status,
      optimizer_converged = model$fit_status$converged,
      active_components_by_stratum = rowSums(model$delta > 1e-8),
      params = iter_params,  # Store all parameters used in this simulation
      N = N,
      simulation_index_within_N = i
    )

    if(!is.null(progress_tick)) {
      progress_tick(sprintf("%s N=%d sim %d/%d", job_label, N, i, n_sims))
    }
  }
  
  return(results)
}

# Main simulation study
