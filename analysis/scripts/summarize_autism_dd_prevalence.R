library(SummarizedExperiment)

# Autism/DD prevalence sweep. Mutational variance depends on prevalence only
# through p/(1-p), so compute each fitted distribution once and rescale exactly.

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
latest <- function(pattern) {
  paths <- sort(list.files(file.path(final_runs_dir, "outputs", "data"),
                           pattern = pattern, full.names = TRUE))
  if (!length(paths)) return(NA_character_)
  tail(paths, 1)
}
autism_file <- get_arg("--autism-model-file", latest("^models_autism_mixsqp_.*\\.Rdata$"))
ddd_file <- get_arg("--ddd-model-file", latest("^models_ddd_mixsqp_.*\\.Rdata$"))
output_dir <- get_arg(
  "--output-dir", file.path(final_runs_dir, "outputs", "derived", "autism_dd")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
for (path in c(autism_file, ddd_file)) if (is.na(path) || !file.exists(path)) stop("Model not found: ", path)
for (source_file in c("io.R", "secondary_analysis_functions.R")) {
  source(if (identical(source_file, "secondary_analysis_functions.R")) file.path(repo_dir, "analysis", "scripts", source_file) else file.path(repo_dir, "R", source_file))
}
autism_env <- new.env(parent = globalenv())
ddd_env <- new.env(parent = globalenv())
load(autism_file, envir = autism_env)
load(ddd_file, envir = ddd_env)

effect_moment <- function(endpoints) {
  out <- numeric(length(endpoints))
  nonzero <- abs(endpoints) > 1e-10
  a <- endpoints[nonzero]
  out[nonzero] <- expm1(2 * a) / (2 * a) - 2 * expm1(a) / a + 1
  out
}
mutvar_constant <- function(model, genetic_data) {
  moments <- effect_moment(model$component_endpoints)
  sum(vapply(seq_len(ncol(model$features)), function(annotation) {
    rows <- model$features[, annotation] == 1
    sum(rows) * 2 * mean(genetic_data$case_rate[rows]) *
      sum(pmax(model$delta[annotation, ], 0) * moments)
  }, numeric(1)))
}
bootstrap_constants <- function(model, genetic_data) {
  vapply(seq_along(model$bootstrap_output$bootstrap_delta), function(iteration) {
    indices <- model$bootstrap_output$bootstrap_samples[, iteration]
    boot_model <- model
    boot_model$features <- model$features[indices, , drop = FALSE]
    boot_model$delta <- model$bootstrap_output$bootstrap_delta[[iteration]]
    mutvar_constant(boot_model, genetic_data[indices, , drop = FALSE])
  }, numeric(1))
}
dataset_constants <- function(data, fitted_models, ptv_name, mis2_name, label) {
  index <- setNames(seq_along(data$loop_vars$names), data$loop_vars$names)
  ptv_index <- index[[ptv_name]]
  mis2_index <- index[[mis2_name]]
  ptv_model <- fitted_models[[1]][[ptv_index]]
  mis2_model <- fitted_models[[1]][[mis2_index]]
  ptv_data <- get_genetic_data(ptv_index, data)$genetic_data
  mis2_data <- get_genetic_data(mis2_index, data)$genetic_data
  # Anchor the displayed curve to the estimates stored with the manuscript
  # models, then apply the exact prevalence scaling. The bootstrap constants
  # below are evaluated analytically to avoid new Monte Carlo noise.
  fitted_prevalence <- data$baseprev * data$prev_factors[1]
  prevalence_scale <- fitted_prevalence / (1 - fitted_prevalence)
  full <- (ptv_model$mutvar_output$total_mutvar * data$ptv_scale_factor +
             mis2_model$mutvar_output$total_mutvar) / prevalence_scale
  ptv_boot <- bootstrap_constants(ptv_model, ptv_data)
  mis2_boot <- bootstrap_constants(mis2_model, mis2_data)
  n_boot <- min(length(ptv_boot), length(mis2_boot))
  list(dataset = label, full = full,
       boot = ptv_boot[seq_len(n_boot)] * data$ptv_scale_factor + mis2_boot[seq_len(n_boot)])
}
autism_constants <- dataset_constants(
  autism_env$autism_data, autism_env$BurdenMLE_DN_models_autism,
  "Combined Probands PTV", "Combined Probands Mis2", "Autism"
)
ddd_constants <- dataset_constants(
  ddd_env$kaplanis_data, ddd_env$BurdenMLE_DN_models_DDD,
  "Proband PTV", "Proband Mis2", "DD"
)
prevalences <- seq(0.005, 0.03, length.out = 20)
make_rows <- function(constants) {
  do.call(rbind, lapply(prevalences, function(prevalence) {
    scale <- prevalence / (1 - prevalence)
    boot <- constants$boot * scale
    data.frame(dataset = constants$dataset, prevalence = prevalence,
               estimate = constants$full * scale,
               lower = quantile(boot, 0.025), upper = quantile(boot, 0.975))
  }))
}
summary <- rbind(make_rows(autism_constants), make_rows(ddd_constants))
write.table(summary, file.path(output_dir, "prevalence_sweep.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Autism/DD prevalence sweep written to", output_dir, "\n")
