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
output_dir <- get_arg(
  "--output-dir", file.path(final_runs_dir, "outputs", "derived", "cohort_ddid")
)
cores <- as.integer(get_arg("--cores", 1))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
for (source_file in c(
  "BurdenMLE_DN.R", "estimate_mutvar.R", "io.R", "likelihoods.R",
  "model.R", "secondary_analysis_functions.R"
)) {
  source(if (identical(source_file, "secondary_analysis_functions.R")) file.path(repo_dir, "analysis", "scripts", source_file) else file.path(repo_dir, "R", source_file))
}
source(file.path(repo_dir, "analysis", "scripts", "model_artifacts.R"))
artifact_info <- load_model_artifact(
  model_manifest, "autism", envir = environment(),
  legacy_path = legacy_model_file
)
model_file <- artifact_info$path
model_index <- setNames(seq_along(autism_data$loop_vars$names), autism_data$loop_vars$names)

exact_mutvar <- function(model, genetic_data, prevalence) {
  estimate_mutvar_trio(model, genetic_data, prevalence)$total_mutvar
}
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
fraction_cases <- function(model, genetic_data, threshold = 5) {
  grid <- posterior_grid(model, genetic_data)
  integrand <- exp(grid$log_rr) * (exp(grid$log_rr) > threshold)
  sum(2 * genetic_data$case_rate * as.vector(grid$probability %*% integrand))
}
bootstrap_model <- function(model, genetic_data, iteration) {
  sample_indices <- model$bootstrap_output$bootstrap_indices[, iteration]
  model_boot <- model
  model_boot$conditional_likelihood <- model$conditional_likelihood[sample_indices, , drop = FALSE]
  model_boot$features <- model$features[sample_indices, , drop = FALSE]
  model_boot$delta <- model$bootstrap_output$bootstrap_delta[[iteration]]
  list(model = model_boot, data = genetic_data[sample_indices, , drop = FALSE])
}

combined_strata <- c("Combined DDID", "Combined Non-DDID")
max_strata <- as.vector(outer(
  c("Male DDID", "Male Non-DDID", "Female DDID", "Female Non-DDID"),
  c("ASC", "SPARK", "GeneDx"), paste
))
strata <- c(combined_strata, max_strata)

summarize_stratum <- function(stratum) {
  is_genedx <- grepl("GeneDx$", stratum)
  model_set_index <- if (is_genedx) 3L else 1L
  prevalence_factor_index <- if (is_genedx) 2L else 1L
  ptv_index <- model_index[[paste(stratum, "PTV")]]
  mis2_index <- model_index[[paste(stratum, "Mis2")]]
  if (is.null(ptv_index) || is.null(mis2_index)) stop("Missing fitted stratum: ", stratum)
  ptv_model <- BurdenMLE_DN_models_autism[[model_set_index]][[ptv_index]]
  mis2_model <- BurdenMLE_DN_models_autism[[model_set_index]][[mis2_index]]
  ptv_data <- get_genetic_data(ptv_index, autism_data)$genetic_data
  mis2_data <- get_genetic_data(mis2_index, autism_data)$genetic_data
  prevalence <- autism_data$loop_vars$prevalences[ptv_index] *
    autism_data$prev_factors[prevalence_factor_index]
  n_boot <- min(length(ptv_model$bootstrap_output$bootstrap_delta),
                length(mis2_model$bootstrap_output$bootstrap_delta))
  mutvar <- exact_mutvar(ptv_model, ptv_data, prevalence) * autism_data$ptv_scale_factor +
    exact_mutvar(mis2_model, mis2_data, prevalence)
  frac <- fraction_cases(ptv_model, ptv_data) * autism_data$ptv_scale_factor +
    fraction_cases(mis2_model, mis2_data)
  mutvar_boot <- frac_boot <- numeric(n_boot)
  for (iteration in seq_len(n_boot)) {
    ptv_boot <- bootstrap_model(ptv_model, ptv_data, iteration)
    mis2_boot <- bootstrap_model(mis2_model, mis2_data, iteration)
    mutvar_boot[iteration] <- exact_mutvar(ptv_boot$model, ptv_boot$data, prevalence) *
      autism_data$ptv_scale_factor + exact_mutvar(mis2_boot$model, mis2_boot$data, prevalence)
    frac_boot[iteration] <- fraction_cases(ptv_boot$model, ptv_boot$data) *
      autism_data$ptv_scale_factor + fraction_cases(mis2_boot$model, mis2_boot$data)
  }
  data.frame(
    stratum = stratum, prevalence = prevalence,
    mutvar_combined = mutvar,
    mutvar_combined_lower = unname(quantile(mutvar_boot, 0.025)),
    mutvar_combined_upper = unname(quantile(mutvar_boot, 0.975)),
    fraccase_RR5 = frac,
    fraccase_RR5_lower = unname(quantile(frac_boot, 0.025)),
    fraccase_RR5_upper = unname(quantile(frac_boot, 0.975)),
    PTV_peneff = ptv_model$penetrance$effective_penetrance,
    PTV_peneff_lower = ptv_model$penetrance$effective_penetrance_CI[1],
    PTV_peneff_upper = ptv_model$penetrance$effective_penetrance_CI[2],
    PTV_effRR = ptv_model$penetrance$effective_penetrance / prevalence,
    PTV_effRR_lower = ptv_model$penetrance$effective_penetrance_CI[1] / prevalence,
    PTV_effRR_upper = ptv_model$penetrance$effective_penetrance_CI[2] / prevalence
  )
}

cat("Summarizing", length(strata), "independently fitted DDID strata...\n")
if (cores > 1 && .Platform$OS.type != "windows") {
  results <- parallel::mclapply(strata, summarize_stratum, mc.cores = cores)
} else {
  results <- lapply(strata, summarize_stratum)
}
ddid_summary <- do.call(rbind, results)
ddid_summary$DDID <- ifelse(grepl("Non-DDID", ddid_summary$stratum), "Non-DDID", "DDID")
ddid_summary$Sex <- sub(" .*", "", ddid_summary$stratum)
ddid_summary$Dataset <- ifelse(
  grepl("(ASC|SPARK|GeneDx)$", ddid_summary$stratum),
  sub(".* ", "", ddid_summary$stratum), "Combined"
)

saveRDS(ddid_summary, file.path(output_dir, "ddid_summary.rds"))
write.table(ddid_summary, file.path(output_dir, "ddid_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("DDID summaries written to", output_dir, "\n")
print(ddid_summary[, c("stratum", "mutvar_combined", "fraccase_RR5", "PTV_effRR")],
      row.names = FALSE)
