library(SummarizedExperiment)

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
repo_dir <- dirname(dirname(dirname(script_file)))
final_runs_dir <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT", unset = file.path(repo_dir, "analysis")), mustWork = FALSE)
model_manifest <- get_arg("--model-manifest", NA_character_)
legacy_model_file <- get_arg("--legacy-model-file", NA_character_)
autism_summary_file <- get_arg(
  "--autism-summary-file",
  file.path(final_runs_dir, "outputs", "derived", "main_autism", "figure2_summary.rds")
)
output_dir <- get_arg(
  "--output-dir", file.path(final_runs_dir, "outputs", "derived", "autism_dd")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(autism_summary_file)) stop("Autism summary not found.")
for (source_file in c(
  "BurdenMLE_DN.R", "estimate_mutvar.R", "io.R", "likelihoods.R",
  "model.R", "secondary_analysis_functions.R"
)) source(if (identical(source_file, "secondary_analysis_functions.R")) file.path(repo_dir, "analysis", "scripts", source_file) else file.path(repo_dir, "R", source_file))
source(file.path(repo_dir, "analysis", "scripts", "model_artifacts.R"))
artifact_info <- load_model_artifact(
  model_manifest, "ddd", envir = environment(),
  legacy_path = legacy_model_file
)
ddd_model_file <- artifact_info$path
autism_summary <- readRDS(autism_summary_file)

ddd_index <- setNames(seq_along(kaplanis_data$loop_vars$names), kaplanis_data$loop_vars$names)
ddd_models <- BurdenMLE_DN_models_DDD[[1]]
get_ddd_model <- function(name) ddd_models[[ddd_index[[name]]]]
get_ddd_data <- function(name) get_genetic_data(ddd_index[[name]], kaplanis_data)$genetic_data
thresholds <- 2:20
prevalence <- kaplanis_data$baseprev * kaplanis_data$prev_factors[1]

posterior_grid <- function(model, genetic_data) {
  component_posterior <- (model$features %*% model$delta) * model$conditional_likelihood
  component_posterior <- component_posterior / rowSums(component_posterior)
  mu_grid <- seq(0.05, 1, by = 1 / model$grid_size)
  probability <- matrix(0, nrow(genetic_data), length(model$component_endpoints) * length(mu_grid))
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
posterior_threshold_summary <- function(model, genetic_data) {
  grid <- posterior_grid(model, genetic_data)
  rr <- exp(grid$log_rr)
  indicators <- vapply(thresholds, function(x) as.numeric(rr > x), numeric(length(rr)))
  list(
    fraction = colSums(2 * genetic_data$case_rate *
                         (grid$probability %*% (rr * indicators))),
    number_genes = colSums(grid$probability %*% indicators)
  )
}
bootstrap_model <- function(model, genetic_data, iteration) {
  sample_indices <- model$bootstrap_output$bootstrap_indices[, iteration]
  model_boot <- model
  model_boot$conditional_likelihood <- model$conditional_likelihood[sample_indices, , drop = FALSE]
  model_boot$features <- model$features[sample_indices, , drop = FALSE]
  model_boot$delta <- model$bootstrap_output$bootstrap_delta[[iteration]]
  list(model = model_boot, data = genetic_data[sample_indices, , drop = FALSE])
}

ptv_model <- get_ddd_model("Proband PTV")
mis2_model <- get_ddd_model("Proband Mis2")
ptv_data <- get_ddd_data("Proband PTV")
mis2_data <- get_ddd_data("Proband Mis2")
ptv_summary <- posterior_threshold_summary(ptv_model, ptv_data)
mis2_summary <- posterior_threshold_summary(mis2_model, mis2_data)
n_boot <- min(length(ptv_model$bootstrap_output$bootstrap_delta),
              length(mis2_model$bootstrap_output$bootstrap_delta))
fraction_boot <- number_boot <- matrix(NA_real_, length(thresholds), n_boot)
fraction_ptv_boot <- fraction_mis2_boot <- matrix(
  NA_real_, length(thresholds), n_boot
)
progress <- txtProgressBar(min = 0, max = n_boot, style = 3)
for (iteration in seq_len(n_boot)) {
  ptv_boot <- bootstrap_model(ptv_model, ptv_data, iteration)
  mis2_boot <- bootstrap_model(mis2_model, mis2_data, iteration)
  ptv_iter <- posterior_threshold_summary(ptv_boot$model, ptv_boot$data)
  mis2_iter <- posterior_threshold_summary(mis2_boot$model, mis2_boot$data)
  fraction_ptv_boot[, iteration] <-
    ptv_iter$fraction * kaplanis_data$ptv_scale_factor
  fraction_mis2_boot[, iteration] <- mis2_iter$fraction
  fraction_boot[, iteration] <-
    fraction_ptv_boot[, iteration] + fraction_mis2_boot[, iteration]
  number_boot[, iteration] <- ptv_iter$number_genes
  setTxtProgressBar(progress, iteration)
}
close(progress)

ddd_fraction <- data.frame(
  threshold = thresholds,
  ptv = ptv_summary$fraction * kaplanis_data$ptv_scale_factor,
  ptv_lower = apply(fraction_ptv_boot, 1, quantile, 0.025),
  ptv_upper = apply(fraction_ptv_boot, 1, quantile, 0.975),
  mis2 = mis2_summary$fraction,
  mis2_lower = apply(fraction_mis2_boot, 1, quantile, 0.025),
  mis2_upper = apply(fraction_mis2_boot, 1, quantile, 0.975),
  estimate = ptv_summary$fraction * kaplanis_data$ptv_scale_factor + mis2_summary$fraction,
  lower = apply(fraction_boot, 1, quantile, 0.025),
  upper = apply(fraction_boot, 1, quantile, 0.975)
)
ddd_number_genes <- data.frame(
  threshold = thresholds, estimate = ptv_summary$number_genes,
  lower = apply(number_boot, 1, quantile, 0.025),
  upper = apply(number_boot, 1, quantile, 0.975)
)

variant_classes <- c("PTV", "Mis2", "Mis1", "Mis0", "Syn")
ddd_mutvar <- do.call(rbind, lapply(variant_classes, function(class) {
  model <- get_ddd_model(paste("Proband", class))
  scale <- if (class == "PTV") kaplanis_data$ptv_scale_factor else 1
  data.frame(
    variant_class = class, estimate = model$mutvar_output$total_mutvar * scale,
    lower = model$mutvar_output$mutvar_CI[1] * scale,
    upper = model$mutvar_output$mutvar_CI[2] * scale
  )
}))
penetrance_classes <- c("PTV", "Mis2", "Mis1")
ddd_effective_rr <- do.call(rbind, lapply(penetrance_classes, function(class) {
  model <- get_ddd_model(paste("Proband", class))
  data.frame(
    variant_class = class,
    estimate = model$penetrance$effective_penetrance / prevalence,
    lower = model$penetrance$effective_penetrance_CI[1] / prevalence,
    upper = model$penetrance$effective_penetrance_CI[2] / prevalence
  )
}))
autism_prevalence <- autism_summary$metadata$prevalence
autism_effective_rr <- transform(
  autism_summary$effective_penetrance,
  estimate = estimate / autism_prevalence,
  lower = lower / autism_prevalence,
  upper = upper / autism_prevalence
)
autism_fraction <- transform(
  autism_summary$fraction_cases,
  estimate = combined, lower = combined_lower, upper = combined_upper
)[c("threshold", "estimate", "lower", "upper")]
ddd_fraction_combined <- ddd_fraction[c("threshold", "estimate", "lower", "upper")]
autism_number_genes <- autism_summary$ptv_number_genes

summary_output <- list(
  metadata = list(ddd_model_file = normalizePath(ddd_model_file),
                  ddd_prevalence = prevalence, autism_prevalence = autism_prevalence),
  mutational_variance = rbind(
    transform(autism_summary$mutational_variance, diagnosis = "Autism"),
    transform(ddd_mutvar, diagnosis = "DD")
  ),
  fraction_cases = rbind(transform(autism_fraction, diagnosis = "Autism"),
                         transform(ddd_fraction_combined, diagnosis = "DD")),
  effective_rr = rbind(transform(autism_effective_rr, diagnosis = "Autism"),
                       transform(ddd_effective_rr, diagnosis = "DD")),
  number_genes = rbind(transform(autism_number_genes, diagnosis = "Autism"),
                       transform(ddd_number_genes, diagnosis = "DD")),
  ddd_fraction_by_class = ddd_fraction
)
saveRDS(summary_output, file.path(output_dir, "figure4_summary.rds"))
for (name in names(summary_output)[-1]) {
  write.table(summary_output[[name]], file.path(output_dir, paste0(name, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
cat("Autism/DD Figure 4 summaries written to", output_dir, "\n")
print(summary_output$mutational_variance, row.names = FALSE)
