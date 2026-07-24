library(SummarizedExperiment)

# Derive the cohort-by-sex quantities used in Figure 2G-H. The estimands and
# bootstrap resampling match visualization_trio.R; only repeated Monte Carlo
# integration is replaced with its exact uniform-mixture expectation.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit)) return(sub(paste0("^", flag, "="), "", hit[1]))
  i <- match(flag, args)
  if (!is.na(i)) return(args[i + 1])
  default
}
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_file <- normalizePath(script_file, mustWork = TRUE)
repo_dir <- dirname(dirname(script_file))
final_runs_dir <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT", unset = file.path(repo_dir, "analysis")), mustWork = FALSE)
model_candidates <- sort(list.files(
  file.path(final_runs_dir, "outputs", "data"),
  pattern = "^models_autism_mixsqp_.*\\.Rdata$", full.names = TRUE
))
model_file <- get_arg(
  "--model-file", if (length(model_candidates)) tail(model_candidates, 1) else NA_character_
)
output_dir <- get_arg(
  "--output-dir", file.path(final_runs_dir, "outputs", "derived", "cohort_ddid")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (is.na(model_file) || !file.exists(model_file)) stop("Main autism model not found.")

for (source_file in c(
  "BurdenMLE_DN.R", "estimate_mutvar.R", "io.R", "likelihoods.R",
  "model.R", "secondary_analysis_functions.R"
)) source(if (identical(source_file, "secondary_analysis_functions.R")) file.path(repo_dir, "analysis", "scripts", source_file) else file.path(repo_dir, "R", source_file))
load(model_file)

model_index <- setNames(seq_along(autism_data$loop_vars$names), autism_data$loop_vars$names)
models <- BurdenMLE_DN_models_autism[[1]]
get_model <- function(name) models[[model_index[[name]]]]
get_data <- function(name) get_genetic_data(model_index[[name]], autism_data)$genetic_data

ptv_name <- "Combined Probands PTV"
mis2_name <- "Combined Probands Mis2"
ptv_model <- get_model(ptv_name)
mis2_model <- get_model(mis2_name)
ptv_data <- get_data(ptv_name)
mis2_data <- get_data(mis2_name)

uniform_exp_mean <- function(endpoint, multiplier) {
  z <- multiplier * endpoint
  ifelse(abs(z) < 1e-12, 1, expm1(z) / z)
}

exact_mutvar <- function(model, genetic_data, prevalence, gamma_scaling_factor = 1) {
  component_moment <-
    uniform_exp_mean(model$component_endpoints, 2 * gamma_scaling_factor) -
    2 * uniform_exp_mean(model$component_endpoints, gamma_scaling_factor) + 1
  annotation_mutvar <- vapply(seq_len(ncol(model$features)), function(annotation) {
    rows <- model$features[, annotation] == 1
    weights <- model$delta[annotation, ]
    weights[weights < 0] <- 0
    weights <- weights / sum(weights)
    sum(rows) * prevalence * mean(2 * genetic_data$case_rate[rows]) *
      sum(weights * component_moment) / (1 - prevalence)
  }, numeric(1))
  sum(annotation_mutvar)
}

posterior_grid <- function(model, genetic_data) {
  component_posterior <- (model$features %*% model$delta) * model$conditional_likelihood
  component_posterior <- component_posterior / rowSums(component_posterior)
  mu_grid <- seq(0.05, 1, by = 1 / model$grid_size)
  probability <- matrix(
    0, nrow(genetic_data), length(model$component_endpoints) * length(mu_grid)
  )
  log_rr <- numeric(ncol(probability))
  counts <- matrix(genetic_data$case_count, nrow(genetic_data), length(mu_grid))
  for (component in seq_along(model$component_endpoints)) {
    columns <- ((component - 1) * length(mu_grid) + 1):(component * length(mu_grid))
    effects <- mu_grid * model$component_endpoints[component]
    likelihood <- dpois(counts, outer(genetic_data$expected_count, exp(effects)))
    likelihood <- likelihood / rowSums(likelihood)
    probability[, columns] <- likelihood * component_posterior[, component]
    log_rr[columns] <- effects
  }
  list(probability = probability, log_rr = log_rr)
}

fraction_cases <- function(model, genetic_data, gamma_scaling_factor, threshold = 5) {
  grid <- posterior_grid(model, genetic_data)
  integrand <- exp(gamma_scaling_factor * grid$log_rr) *
    (exp(grid$log_rr) > threshold)
  sum(2 * genetic_data$case_rate * as.vector(grid$probability %*% integrand))
}

bootstrap_model <- function(model, genetic_data, iteration) {
  sample_indices <- model$bootstrap_output$bootstrap_samples[, iteration]
  model_boot <- model
  model_boot$conditional_likelihood <-
    model$conditional_likelihood[sample_indices, , drop = FALSE]
  model_boot$features <- model$features[sample_indices, , drop = FALSE]
  model_boot$delta <- model$bootstrap_output$bootstrap_delta[[iteration]]
  list(model = model_boot, data = genetic_data[sample_indices, , drop = FALSE])
}

if (!identical(
  ptv_model$bootstrap_output$bootstrap_samples,
  mis2_model$bootstrap_output$bootstrap_samples
)) stop("PTV and Mis2 bootstrap samples must be paired.")
n_boot <- min(
  length(ptv_model$bootstrap_output$bootstrap_delta),
  length(mis2_model$bootstrap_output$bootstrap_delta)
)

scenario_rows <- expand.grid(
  study = c("SPARK", "ASC", "GeneDx"),
  sex = c("Male", "Female", "Both"),
  stringsAsFactors = FALSE
)
scenario_rows$model_name <- ifelse(
  scenario_rows$sex == "Both",
  paste("Combined", scenario_rows$study, "PTV"),
  paste(scenario_rows$sex, scenario_rows$study, "PTV")
)
if (anyNA(model_index[scenario_rows$model_name])) {
  stop("Missing cohort/sex entries: ",
       paste(scenario_rows$model_name[is.na(model_index[scenario_rows$model_name])], collapse = ", "))
}

cat("Estimating cohort/sex scaling factors...\n")
scenario_rows$gamma_scaling_factor <- vapply(scenario_rows$model_name, function(name) {
  binomial_analysis(
    model = ptv_model, genetic_data_total = ptv_data,
    genetic_data_subsample = get_data(name)
  )$MLE
}, numeric(1))
scenario_rows$prevalence <- vapply(seq_len(nrow(scenario_rows)), function(i) {
  factor_index <- if (scenario_rows$study[i] %in% c("SPARK", "ASC")) 1 else 2
  autism_data$loop_vars$prevalences[model_index[[scenario_rows$model_name[i]]]] *
    autism_data$prev_factors[factor_index]
}, numeric(1))

results <- vector("list", nrow(scenario_rows))
progress <- txtProgressBar(min = 0, max = nrow(scenario_rows), style = 3)
for (scenario_index in seq_len(nrow(scenario_rows))) {
  scaling <- scenario_rows$gamma_scaling_factor[scenario_index]
  prevalence <- scenario_rows$prevalence[scenario_index]
  mutvar_ptv <- exact_mutvar(ptv_model, ptv_data, prevalence, scaling) *
    autism_data$ptv_scale_factor
  mutvar_mis2 <- exact_mutvar(mis2_model, mis2_data, prevalence, scaling)
  frac_ptv <- fraction_cases(ptv_model, ptv_data, scaling) * autism_data$ptv_scale_factor
  frac_mis2 <- fraction_cases(mis2_model, mis2_data, scaling)

  mutvar_boot <- numeric(n_boot)
  frac_boot <- numeric(n_boot)
  for (iteration in seq_len(n_boot)) {
    ptv_boot <- bootstrap_model(ptv_model, ptv_data, iteration)
    mis2_boot <- bootstrap_model(mis2_model, mis2_data, iteration)
    mutvar_boot[iteration] <-
      exact_mutvar(ptv_boot$model, ptv_boot$data, prevalence, scaling) *
        autism_data$ptv_scale_factor +
      exact_mutvar(mis2_boot$model, mis2_boot$data, prevalence, scaling)
    frac_boot[iteration] <-
      fraction_cases(ptv_boot$model, ptv_boot$data, scaling) *
        autism_data$ptv_scale_factor +
      fraction_cases(mis2_boot$model, mis2_boot$data, scaling)
  }
  results[[scenario_index]] <- cbind(
    scenario_rows[scenario_index, ],
    data.frame(
      mutvar_ptv = mutvar_ptv, mutvar_mis2 = mutvar_mis2,
      mutvar_combined = mutvar_ptv + mutvar_mis2,
      mutvar_combined_lower = unname(quantile(mutvar_boot, 0.025)),
      mutvar_combined_upper = unname(quantile(mutvar_boot, 0.975)),
      fraccase_RR5 = frac_ptv + frac_mis2,
      fraccase_RR5_lower = unname(quantile(frac_boot, 0.025)),
      fraccase_RR5_upper = unname(quantile(frac_boot, 0.975))
    )
  )
  setTxtProgressBar(progress, scenario_index)
}
close(progress)
cohort_sex_summary <- do.call(rbind, results)
cohort_sex_summary$study[cohort_sex_summary$study == "GeneDx"] <- "GeneDX"

saveRDS(cohort_sex_summary, file.path(output_dir, "cohort_sex_summary.rds"))
write.table(
  cohort_sex_summary, file.path(output_dir, "cohort_sex_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
cat("Cohort/sex summary written to", output_dir, "\n")
print(cohort_sex_summary[, c(
  "study", "sex", "prevalence", "gamma_scaling_factor",
  "mutvar_combined", "fraccase_RR5"
)], row.names = FALSE)
