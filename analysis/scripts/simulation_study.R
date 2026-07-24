library(ggplot2)
library(data.table)
library(patchwork)
library(bayess)

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default) {
  equals_match <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(equals_match) > 0) {
    return(sub(paste0("^", flag, "="), "", equals_match[1]))
  }
  flag_index <- match(flag, args)
  if (!is.na(flag_index)) {
    if (flag_index == length(args) || startsWith(args[flag_index + 1], "--")) {
      stop("Missing value after ", flag)
    }
    return(args[flag_index + 1])
  }
  default
}

if ("--help" %in% args) {
  cat(
    "Options:\n",
    "  --optimizer <EM|mixsqp>  Model optimizer (default: EM)\n",
    "  --run-date <label>        Output date label\n",
    "  --seed <integer>          Simulation seed (default: 20260715)\n",
    sep = ""
  )
  quit(status = 0)
}

optimizer <- get_arg("--optimizer", "EM")
if (tolower(optimizer) == "em") {
  optimizer <- "EM"
} else if (tolower(optimizer) == "mixsqp") {
  optimizer <- "mixsqp"
} else {
  stop("--optimizer must be EM or mixsqp.")
}
optimizer_label <- tolower(optimizer)

script_file <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
)
script_file <- normalizePath(script_file, mustWork = TRUE)
repo_dir <- dirname(dirname(script_file))
final_runs_dir <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT", unset = file.path(repo_dir, "analysis")), mustWork = FALSE)
output_data_dir <- file.path(final_runs_dir, "outputs", "data")

dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)

run_date <- get_arg("--run-date", format(Sys.Date(), "%b%d_%y"))
if (!grepl("^[A-Za-z0-9_-]+$", run_date)) {
  stop("--run-date may contain only letters, numbers, underscores, and hyphens.")
}

seed_str <- get_arg("--seed", Sys.getenv("SIM_RNG_SEED"))
if (seed_str == "") seed_str <- "20260715"
rng_seed <- as.integer(seed_str)
if (is.na(rng_seed)) stop("--seed must be an integer.")
set.seed(rng_seed)

effect_size_samples_file <- file.path(output_data_dir,
  sprintf("effect_size_samples_%s_seed%s.RData", run_date, rng_seed))
simulation_results_file <- file.path(output_data_dir,
  sprintf("simulation_results_%s_%s_seed%s.RData",
          optimizer_label, run_date, rng_seed))

source_files <- c(
  "BurdenMLE_DN.R", "EM.R", "MixSQP.R", "estimate_mutvar.R",
  "likelihoods.R", "model.R", "io.R", "secondary_analysis_functions.R"
)
for (source_file in source_files) {
  source(if (identical(source_file, "secondary_analysis_functions.R")) file.path(repo_dir, "analysis", "scripts", source_file) else file.path(repo_dir, "R", source_file))
}
source(file.path(repo_dir, "analysis", "scripts", "set_up_asc.R"))
# Function to generate effect sizes from different distributions
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
  # Calculate numerator: mean((exp(x)-1)*exp(x))
  peneff_numerator <- mean((exp(true_effects) - 1) * exp(true_effects))
  
  # Calculate denominator: mean(exp(x)-1)
  peneff_denominator <- mean(exp(true_effects) - 1)
  
  # Calculate effective penetrance
  peneff <- prevalence * (peneff_numerator/peneff_denominator)
  return(peneff)
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

    final_log_likelihood <- if (optimizer == "EM") {
      tail(stats::na.omit(model$ll), 1)
    } else {
      model$ll
    }
    optimizer_status <- if (optimizer == "EM") {
      if (sum(!is.na(model$ll)) < 10000) "tolerance reached" else "maximum iterations reached"
    } else {
      paste(vapply(model$mixsqp_output, function(x) x$status, character(1)),
            collapse = "; ")
    }

    #get the estimated fraction of cases with large effect mutations, using the same threshold
    est_fraccases <- get_fraccase(model, sim_data, iter_params$RR_threshold)

    # Store estimates and compact diagnostics, but not the full fitted model.
    results[[i]] <- list(
      true_effects = true_effects,
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
      optimizer_iterations = if (optimizer == "EM") sum(!is.na(model$ll)) else NA_integer_,
      optimizer_status = optimizer_status,
      optimizer_converged = if (optimizer == "EM") {
        sum(!is.na(model$ll)) < 10000
      } else {
        all(vapply(model$mixsqp_output,
                   function(x) x$status == "converged to optimal solution",
                   logical(1)))
      },
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
run_simulation_study <- function() {
  # Get empirical data structure from autism analysis
  processed_input <- get_genetic_data(5, autism_data)  # Using PTV data as reference
  
  # Sample sizes to simulate
  sample_sizes <- c(1000, 5000, 15000, 38680, 50000)
  
  # Number of bins for half_uniform
  n_bins <- 10  # Back to 10 components
  
  # Define simulation parameters
  sim_params <- list(
    RR_threshold = 5,
    prevalence = prev_weighted,  # Use weighted prevalence across cohorts
    
    # Distribution-specific parameters
    null = list(),
    
    half_uniform = list(
      endpoints = seq(0, log(1/0.01), length.out = n_bins),
      weights = function() bayess::rdirichlet(1, rep(1, n_bins))[1,]
    ),
    
    half_infinitesimal = list(
      # Sample sd uniformly between 0.5 and 1 for each simulation
      sd = function() runif(1, 0.5, 1)
    ),
    
    point_normal = list(
      # Sample pi0 between 0.85 and 0.95
      pi0 = function() runif(1, 0.85, 0.95),
      # Sample mean between 1.5 and 2.5 to stay well above 0
      mean = function() runif(1, 1.5, 2.5),
      # Sample sd between 0.2 and 0.4 to minimize negative values
      sd = function() runif(1, 0.2, 0.4)
    ),
    
    oligogenic = list(
      # Sample prop_causal between 0.001 and 0.01 (0.1% to 1% of genes)
      prop_causal = function() runif(1, 0.001, 0.01),
      # Sample effect size uniformly between log(5) and log(50)
      effect_size = function() runif(1, log(5), log(50))
    )
  )
  
  # Run simulations for each distribution and sample size
  results <- list()
  distributions <- c("null", "half_uniform", "half_infinitesimal", "point_normal", "oligogenic")
  n_sims_per_job <- 100

  total_sims <- length(distributions) * length(sample_sizes) * n_sims_per_job
  use_progress_pkg <- requireNamespace("progress", quietly = TRUE)

  if(use_progress_pkg) {
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent eta: :eta | :current/:total | :message",
      total = total_sims,
      clear = FALSE,
      width = 90
    )
    progress_tick <- function(message) {
      pb$tick(tokens = list(message = message))
    }
  } else {
    pb <- txtProgressBar(min = 0, max = total_sims, style = 3)
    progress_counter <- 0
    progress_tick <- function(message) {
      progress_counter <<- progress_counter + 1
      setTxtProgressBar(pb, progress_counter)
    }
    cat("progress package not installed; using base R progress bar.\n")
  }
  
  # Generate effect size samples for visualization using number of genes from processed data
  n_samples <- nrow(processed_input$genetic_data)  # Use actual number of genes from autism data
  effect_size_samples <- list()
  
  # Sample from each distribution once with random parameters
  for(dist_name in distributions) {
    if(dist_name == "null") {
      samples <- rep(0, n_samples)
    } else {
      samples <- switch(dist_name,
        "half_uniform" = {
          weights <- sim_params[[dist_name]]$weights()
          component <- sample(1:length(sim_params[[dist_name]]$endpoints), n_samples, 
                            prob = weights, replace = TRUE)
          sapply(1:n_samples, function(i) 
            runif(1, min=0, max=sim_params[[dist_name]]$endpoints[component[i]]))
        },
        "half_infinitesimal" = {
          sd <- sim_params[[dist_name]]$sd()
          abs(rnorm(n_samples, mean=0, sd=sd))
        },
        "point_normal" = {
          pi0 <- sim_params[[dist_name]]$pi0()
          mean <- sim_params[[dist_name]]$mean()
          sd <- sim_params[[dist_name]]$sd()
          is_zero <- rbinom(n_samples, 1, pi0)
          effects <- abs(rnorm(n_samples, mean=mean, sd=sd))
          effects[is_zero == 1] <- 0
          effects
        },
        "oligogenic" = {
          # Majority zero, few large effects
          effects <- rep(0, n_samples)
          n_causal <- round(n_samples * sim_params[[dist_name]]$prop_causal())
          causal_idx <- sample(1:n_samples, n_causal)
          effects[causal_idx] <- sim_params[[dist_name]]$effect_size()  # Single effect size for all causal genes
          effects
        }
      )
    }
    effect_size_samples[[dist_name]] <- samples
  }
  
  # Convert to data frame
  effect_size_df <- data.frame(
    effect_size = unlist(effect_size_samples),
    distribution = rep(names(effect_size_samples), each = n_samples)
  )
  
  # Save effect size samples for visualization
  save(effect_size_df, file = effect_size_samples_file)
  
  for(dist in distributions) {
    results[[dist]] <- list()  # Initialize list for this distribution
    sim_count <- 1  # Counter for simulations within this distribution
    
    for(N in sample_sizes) {
      cat(sprintf("\nRunning simulations for %s distribution with N=%d...\n", dist, N))
      
      # Create a copy of the base parameters
      base_params <- sim_params[[dist]]
      
      # Run simulations for this N
      N_results <- run_simulation(
        distribution = dist,
        params = c(base_params, list(prevalence = sim_params$prevalence, RR_threshold = sim_params$RR_threshold)),
        empirical_data = processed_input$genetic_data,
        empirical_features = processed_input$features,
        N = N,
        n_sims = n_sims_per_job,
        optimizer = optimizer,
        progress_tick = progress_tick,
        job_label = dist
      )
      
      # Add results to the list
      results[[dist]] <- c(results[[dist]], N_results)
    }
  }

  if(use_progress_pkg) {
    cat("\n")
  } else {
    close(pb)
    cat("\n")
  }
  
  # Save results
  simulation_metadata <- list(
    optimizer = optimizer,
    rng_seed = rng_seed,
    run_date = run_date,
    sample_sizes = sample_sizes,
    distributions = distributions,
    simulations_per_distribution_and_N = n_sims_per_job,
    component_endpoints = seq(0, log(1 / 0.01), length.out = n_bins),
    optimizer_rng_isolated = TRUE
  )
  save(results, simulation_metadata,
       file = simulation_results_file)

  cat(sprintf("Saved effect size samples: %s\n", effect_size_samples_file))
  cat(sprintf("Saved simulation results: %s\n", simulation_results_file))
  
  return(results)
}

# Run the simulation study
results <- run_simulation_study()
