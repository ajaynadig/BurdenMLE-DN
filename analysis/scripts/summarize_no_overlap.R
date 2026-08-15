library(SummarizedExperiment)

# Autism/DD sensitivity after removing GeneDx samples overlapping the DD study.

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
data_dir <- file.path(final_runs_dir, "outputs", "data")
latest_no_overlap <- function() {
  candidates <- list.files(
    data_dir,
    pattern = "^models_autism_no_overlap_mixsqp_.*\\.Rdata$",
    full.names = TRUE
  )
  if (!length(candidates)) {
    stop("No dated no-overlap MixSQP model file was found in ", data_dir)
  }
  candidates[which.max(file.info(candidates)$mtime)]
}
model_file <- get_arg("--model-file", NA_character_)
if (is.na(model_file)) model_file <- latest_no_overlap()
figure4_summary_file <- get_arg(
  "--figure4-summary-file", file.path(final_runs_dir, "outputs", "derived", "autism_dd", "figure4_summary.rds")
)
output_dir <- get_arg(
  "--output-dir", file.path(final_runs_dir, "outputs", "derived", "sensitivities")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
for (path in c(model_file, figure4_summary_file)) if (!file.exists(path)) stop("Required input not found: ", path)
for (source_file in c("io.R", "secondary_analysis_functions.R")) source(if (identical(source_file, "secondary_analysis_functions.R")) file.path(repo_dir, "analysis", "scripts", source_file) else file.path(repo_dir, "R", source_file))
load(model_file)
main_comparison <- readRDS(figure4_summary_file)
data <- autism_data_NoKaplanis
fitted_models <- BurdenMLE_DN_models_autism_NoKaplanis
index <- setNames(seq_along(data$loop_vars$names), data$loop_vars$names)
models <- fitted_models[[1]]
get_model <- function(name) models[[index[[name]]]]
get_data <- function(name) get_genetic_data(index[[name]], data)$genetic_data
thresholds <- 2:20

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
threshold_summary <- function(model, genetic_data) {
  grid <- posterior_grid(model, genetic_data)
  rr <- exp(grid$log_rr)
  indicators <- vapply(thresholds, function(x) as.numeric(rr > x), numeric(length(rr)))
  list(
    fraction = colSums(2 * genetic_data$case_rate * (grid$probability %*% (rr * indicators))),
    genes = colSums(grid$probability %*% indicators)
  )
}
bootstrap_model <- function(model, genetic_data, iteration) {
  indices <- model$bootstrap_output$bootstrap_indices[, iteration]
  boot <- model
  boot$conditional_likelihood <- model$conditional_likelihood[indices, , drop = FALSE]
  boot$features <- model$features[indices, , drop = FALSE]
  boot$delta <- model$bootstrap_output$bootstrap_delta[[iteration]]
  list(model = boot, data = genetic_data[indices, , drop = FALSE])
}

ptv_model <- get_model("Combined Probands PTV")
mis2_model <- get_model("Combined Probands Mis2")
ptv_data <- get_data("Combined Probands PTV")
mis2_data <- get_data("Combined Probands Mis2")
ptv <- threshold_summary(ptv_model, ptv_data)
mis2 <- threshold_summary(mis2_model, mis2_data)
n_boot <- min(length(ptv_model$bootstrap_output$bootstrap_delta),
              length(mis2_model$bootstrap_output$bootstrap_delta))
fraction_boot <- gene_boot <- matrix(NA_real_, length(thresholds), n_boot)
progress <- txtProgressBar(min = 0, max = n_boot, style = 3)
for (iteration in seq_len(n_boot)) {
  ptv_boot <- bootstrap_model(ptv_model, ptv_data, iteration)
  mis2_boot <- bootstrap_model(mis2_model, mis2_data, iteration)
  ptv_iter <- threshold_summary(ptv_boot$model, ptv_boot$data)
  mis2_iter <- threshold_summary(mis2_boot$model, mis2_boot$data)
  fraction_boot[, iteration] <- ptv_iter$fraction * data$ptv_scale_factor + mis2_iter$fraction
  gene_boot[, iteration] <- ptv_iter$genes
  setTxtProgressBar(progress, iteration)
}
close(progress)

autism_fraction <- data.frame(
  threshold = thresholds,
  estimate = ptv$fraction * data$ptv_scale_factor + mis2$fraction,
  lower = apply(fraction_boot, 1, quantile, 0.025),
  upper = apply(fraction_boot, 1, quantile, 0.975), diagnosis = "Autism"
)
autism_genes <- data.frame(
  threshold = thresholds, estimate = ptv$genes,
  lower = apply(gene_boot, 1, quantile, 0.025),
  upper = apply(gene_boot, 1, quantile, 0.975), diagnosis = "Autism"
)
variant_classes <- c("PTV", "Mis2", "Mis1", "Mis0", "Syn")
autism_mutvar <- do.call(rbind, lapply(variant_classes, function(class) {
  model <- get_model(paste("Combined Probands", class))
  scale <- if (class == "PTV") data$ptv_scale_factor else 1
  data.frame(variant_class = class, estimate = model$mutvar_output$total_mutvar * scale,
             lower = model$mutvar_output$mutvar_CI[1] * scale,
             upper = model$mutvar_output$mutvar_CI[2] * scale, diagnosis = "Autism")
}))
prevalence <- data$baseprev * data$prev_factors[1]
autism_rr <- do.call(rbind, lapply(c("PTV", "Mis2", "Mis1"), function(class) {
  model <- get_model(paste("Combined Probands", class))
  data.frame(variant_class = class,
             estimate = model$penetrance$effective_penetrance / prevalence,
             lower = model$penetrance$effective_penetrance_CI[1] / prevalence,
             upper = model$penetrance$effective_penetrance_CI[2] / prevalence,
             diagnosis = "Autism")
}))

dd_rows <- function(table) table[table$diagnosis == "DD", ]
summary_output <- list(
  metadata = list(autism_model_file = normalizePath(model_file), optimizer = "mixsqp",
                  note = "Dedicated no-overlap autism fit compared with current DD fit"),
  mutational_variance = rbind(autism_mutvar, dd_rows(main_comparison$mutational_variance)),
  fraction_cases = rbind(autism_fraction, dd_rows(main_comparison$fraction_cases)),
  effective_rr = rbind(autism_rr, dd_rows(main_comparison$effective_rr)),
  number_genes = rbind(autism_genes, dd_rows(main_comparison$number_genes))
)
saveRDS(summary_output, file.path(output_dir, "no_overlap_summary.rds"))
for (name in names(summary_output)[-1]) write.table(
  summary_output[[name]], file.path(output_dir, paste0("no_overlap_", name, ".tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
cat("No-overlap sensitivity summary written to", output_dir, "\n")
