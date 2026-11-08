library(ggplot2)
library(data.table)
library(patchwork)


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
repo_dir <- dirname(dirname(dirname(script_file)))
source(file.path(repo_dir, "analysis/scripts/simulation_functions.R"))
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
      weights = function() draw_dirichlet(rep(1, n_bins))
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
