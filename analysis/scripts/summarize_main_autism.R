# Generate the derived numerical results used by the central autism panels of
# Figure 2. This script performs secondary analysis only; plotting is handled
# by make_figure2_main_autism.R.

library(data.table)
library(SummarizedExperiment)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  equals_match <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(equals_match)) {
    return(sub(paste0("^", flag, "="), "", equals_match[1]))
  }
  flag_index <- match(flag, args)
  if (!is.na(flag_index)) return(args[flag_index + 1])
  default
}

script_file <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
)
script_file <- normalizePath(script_file, mustWork = TRUE)
repo_dir <- dirname(dirname(script_file))
final_runs_dir <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT", unset = file.path(repo_dir, "analysis")), mustWork = FALSE)

model_candidates <- sort(list.files(
  file.path(final_runs_dir, "outputs", "data"),
  pattern = "^models_autism_mixsqp_.*\\.Rdata$",
  full.names = TRUE
))
model_file <- get_arg(
  "--model-file",
  if (length(model_candidates)) tail(model_candidates, 1) else NA_character_
)
output_dir <- get_arg(
  "--output-dir",
  file.path(final_runs_dir, "outputs", "derived", "main_autism")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (is.na(model_file) || !file.exists(model_file))
  stop("No main autism MixSQP model file was found: ", model_file)

source_files <- c(
  "BurdenMLE_DN.R", "estimate_mutvar.R", "io.R", "likelihoods.R",
  "model.R", "secondary_analysis_functions.R"
)
for (source_file in source_files) {
  source(if (identical(source_file, "secondary_analysis_functions.R")) file.path(repo_dir, "analysis", "scripts", source_file) else file.path(repo_dir, "R", source_file))
}
load(model_file)

if (!exists("autism_data") || !exists("BurdenMLE_DN_models_autism")) {
  stop("The model file must contain autism_data and BurdenMLE_DN_models_autism.")
}

model_index <- setNames(
  seq_along(autism_data$loop_vars$names),
  autism_data$loop_vars$names
)
required_models <- c(
  "Combined Probands PTV", "Combined Probands Mis2",
  "Combined Probands Mis1", "Combined Probands Mis0",
  "Combined Probands Syn", "Combined SPARK PTV", "Combined ASC PTV",
  "Combined GeneDx PTV"
)
if (anyNA(model_index[required_models])) {
  stop(
    "Required models are missing: ",
    paste(required_models[is.na(model_index[required_models])], collapse = ", ")
  )
}

prevalence_index <- 1L
prevalence <- autism_data$baseprev * autism_data$prev_factors[prevalence_index]
models <- BurdenMLE_DN_models_autism[[prevalence_index]]

get_model <- function(name) models[[model_index[[name]]]]
get_data <- function(name) {
  get_genetic_data(model_index[[name]], autism_data)$genetic_data
}

# Figure 2A describes the component autism cohorts.
summarize_cohort <- function(dataset, model_data) {
  genetic_data <- model_data$genetic_data
  constrained <- model_data$features[, 1] == 1 | model_data$features[, 2] == 1
  observed <- sum(genetic_data$case_count[constrained])
  expected <- sum(genetic_data$expected_count[constrained])
  data.frame(
    dataset = dataset,
    n_proband = genetic_data$N[1],
    n_constrained_ptv = observed,
    expected_constrained_ptv = expected,
    rr_constrained_ptv = observed / expected
  )
}

cohort_names <- c("SPARK", "ASC", "GeneDx")
cohort_model_names <- paste("Combined", cohort_names, "PTV")
cohort_table <- do.call(rbind, lapply(seq_along(cohort_names), function(i) {
  summarize_cohort(
    cohort_names[i],
    get_genetic_data(model_index[[cohort_model_names[i]]], autism_data)
  )
}))

# DDD is stored as a dataset column but is not a fitted standalone loop entry.
# Create the equivalent loop specification and pass it through exactly the same
# get_genetic_data() path as the other cohorts.
ddd_data <- autism_data
ddd_index <- length(ddd_data$loop_vars$names) + 1L
combined_index <- model_index[["Combined Probands PTV"]]
ddd_data$loop_vars$names <- c(ddd_data$loop_vars$names, "Combined DDD PTV")
ddd_data$loop_vars$datasets[[ddd_index]] <-
  colnames(ddd_data$counts) == "DDD ASD"
ddd_data$loop_vars$N_subset[[ddd_index]] <-
  SummarizedExperiment::colData(ddd_data$counts)["DDD ASD", "N_Proband"]
ddd_data$loop_vars$subsets[[ddd_index]] <-
  ddd_data$loop_vars$subsets[[combined_index]]
ddd_data$loop_vars$mutation_rate[[ddd_index]] <-
  ddd_data$loop_vars$mutation_rate[[combined_index]]
ddd_data$loop_vars$prevalences[[ddd_index]] <-
  ddd_data$loop_vars$prevalences[[combined_index]]
ddd_data$loop_vars$fit_model[[ddd_index]] <- FALSE
ddd_row <- summarize_cohort("DDD", get_genetic_data(ddd_index, ddd_data))

all_row <- summarize_cohort(
  "All", get_genetic_data(model_index[["Combined Probands PTV"]], autism_data)
)
cohort_table <- rbind(
  cohort_table,
  ddd_row,
  all_row
)

# Compute the complete discrete approximation to each gene's posterior once.
# Subsequent posterior expectations for all thresholds are matrix products,
# avoiding repeated Poisson likelihood evaluation for every estimand.
posterior_grid <- function(model, genetic_data) {
  component_weights <- model$features %*% model$delta
  component_posteriors <- component_weights * model$conditional_likelihood
  component_posteriors <- component_posteriors / rowSums(component_posteriors)

  mu_grid <- seq(0.05, 1, by = 1 / model$grid_size)
  no_genes <- nrow(genetic_data)
  no_components <- length(model$component_endpoints)
  point_posteriors <- matrix(
    0, nrow = no_genes, ncol = no_components * length(mu_grid)
  )
  log_rr <- numeric(ncol(point_posteriors))
  case_count_matrix <- matrix(
    genetic_data$case_count,
    nrow = no_genes,
    ncol = length(mu_grid)
  )

  for (component in seq_len(no_components)) {
    columns <- ((component - 1) * length(mu_grid) + 1):
      (component * length(mu_grid))
    component_log_rr <- mu_grid * model$component_endpoints[component]
    rate <- outer(genetic_data$expected_count, exp(component_log_rr))
    conditional_points <- dpois(case_count_matrix, rate)
    conditional_points <- conditional_points / rowSums(conditional_points)
    point_posteriors[, columns] <-
      conditional_points * component_posteriors[, component]
    log_rr[columns] <- component_log_rr
  }

  list(probability = point_posteriors, log_rr = log_rr)
}

posterior_summaries <- function(model, genetic_data, thresholds) {
  grid <- posterior_grid(model, genetic_data)
  rate_ratio <- exp(grid$log_rr)
  threshold_functions <- vapply(
    thresholds,
    function(threshold) rate_ratio * (rate_ratio > threshold),
    numeric(length(rate_ratio))
  )
  list(
    mutvar_gene = (
      (prevalence * 2 * genetic_data$case_rate) / (1 - prevalence)
    ) * as.vector(grid$probability %*% (rate_ratio - 1)^2),
    fraction_cases = as.vector(
      colSums((2 * genetic_data$case_rate) *
                (grid$probability %*% threshold_functions))
    ),
    number_genes = as.vector(
      colSums(grid$probability %*%
                vapply(
                  thresholds,
                  function(threshold) as.numeric(rate_ratio > threshold),
                  numeric(length(rate_ratio))
                ))
    )
  )
}

bootstrap_model <- function(model, genetic_data, iteration) {
  sample_indices <- model$bootstrap_output$bootstrap_indices[, iteration]
  model_boot <- model
  model_boot$conditional_likelihood <-
    model$conditional_likelihood[sample_indices, , drop = FALSE]
  model_boot$features <- model$features[sample_indices, , drop = FALSE]
  model_boot$delta <- model$bootstrap_output$bootstrap_delta[[iteration]]
  list(model = model_boot, data = genetic_data[sample_indices, , drop = FALSE])
}

ptv_name <- "Combined Probands PTV"
mis2_name <- "Combined Probands Mis2"
ptv_model <- get_model(ptv_name)
mis2_model <- get_model(mis2_name)
ptv_data <- get_data(ptv_name)
mis2_data <- get_data(mis2_name)
thresholds <- 2:20

if (!identical(
  ptv_model$bootstrap_output$bootstrap_samples,
  mis2_model$bootstrap_output$bootstrap_samples
)) {
  stop("PTV and Mis2 bootstrap samples are not paired.")
}

message("Computing full-data posterior summaries...")
ptv_summary <- posterior_summaries(ptv_model, ptv_data, thresholds)
mis2_summary <- posterior_summaries(mis2_model, mis2_data, thresholds)
figure3_thresholds <- c(2, 5, 10, 20)
ptv_grid_figure3 <- posterior_grid(ptv_model, ptv_data)
ptv_rr_figure3 <- exp(ptv_grid_figure3$log_rr)
ptv_gene_rr_probability <- data.frame(
  gene_id = rownames(ptv_data),
  ptv_grid_figure3$probability %*%
    vapply(
      figure3_thresholds,
      function(threshold) as.numeric(ptv_rr_figure3 > threshold),
      numeric(length(ptv_rr_figure3))
    ),
  check.names = FALSE
)
names(ptv_gene_rr_probability)[-1] <- paste0("rr_greater_", figure3_thresholds)
combined_mutvar_gene <-
  ptv_summary$mutvar_gene * autism_data$ptv_scale_factor +
  mis2_summary$mutvar_gene
names(combined_mutvar_gene) <- rownames(ptv_data)
combined_mutvar_gene <- sort(combined_mutvar_gene, decreasing = TRUE)

n_boot <- min(
  length(ptv_model$bootstrap_output$bootstrap_delta),
  length(mis2_model$bootstrap_output$bootstrap_delta)
)
polygenicity_boot <- matrix(NA_real_, nrow(ptv_data), n_boot)
polygenicity_boot_absolute <- matrix(NA_real_, nrow(ptv_data), n_boot)
fraction_cases_boot <- matrix(NA_real_, nrow = length(thresholds), ncol = n_boot)
fraction_cases_ptv_boot <- matrix(NA_real_, nrow = length(thresholds), ncol = n_boot)
fraction_cases_mis2_boot <- matrix(NA_real_, nrow = length(thresholds), ncol = n_boot)
number_genes_boot <- matrix(NA_real_, nrow = length(thresholds), ncol = n_boot)

message("Computing ", n_boot, " paired bootstrap summaries...")
progress <- txtProgressBar(min = 0, max = n_boot, style = 3)
for (iteration in seq_len(n_boot)) {
  ptv_boot <- bootstrap_model(ptv_model, ptv_data, iteration)
  mis2_boot <- bootstrap_model(mis2_model, mis2_data, iteration)
  ptv_boot_summary <- posterior_summaries(
    ptv_boot$model, ptv_boot$data, thresholds
  )
  mis2_boot_summary <- posterior_summaries(
    mis2_boot$model, mis2_boot$data, thresholds
  )

  combined_boot <-
    ptv_boot_summary$mutvar_gene * autism_data$ptv_scale_factor +
    mis2_boot_summary$mutvar_gene
  combined_boot <- sort(combined_boot, decreasing = TRUE)
  cumulative_boot <- cumsum(combined_boot)
  polygenicity_boot_absolute[, iteration] <- cumulative_boot
  polygenicity_boot[, iteration] <- cumulative_boot / sum(combined_boot)
  fraction_cases_ptv_boot[, iteration] <-
    ptv_boot_summary$fraction_cases * autism_data$ptv_scale_factor
  fraction_cases_mis2_boot[, iteration] <- mis2_boot_summary$fraction_cases
  fraction_cases_boot[, iteration] <-
    fraction_cases_ptv_boot[, iteration] +
    fraction_cases_mis2_boot[, iteration]
  number_genes_boot[, iteration] <- ptv_boot_summary$number_genes
  setTxtProgressBar(progress, iteration)
}
close(progress)

polygenicity <- data.frame(
  number_genes = seq_along(combined_mutvar_gene),
  cumulative_mutvar = cumsum(combined_mutvar_gene) / sum(combined_mutvar_gene),
  lower = apply(polygenicity_boot, 1, quantile, 0.025),
  upper = apply(polygenicity_boot, 1, quantile, 0.975),
  cumulative_mutvar_absolute = cumsum(combined_mutvar_gene),
  lower_absolute = apply(polygenicity_boot_absolute, 1, quantile, 0.025),
  upper_absolute = apply(polygenicity_boot_absolute, 1, quantile, 0.975)
)
half_index <- which(polygenicity$cumulative_mutvar >= 0.5)[1]

fraction_cases <- data.frame(
  threshold = thresholds,
  ptv = ptv_summary$fraction_cases * autism_data$ptv_scale_factor,
  ptv_lower = apply(fraction_cases_ptv_boot, 1, quantile, 0.025),
  ptv_upper = apply(fraction_cases_ptv_boot, 1, quantile, 0.975),
  mis2 = mis2_summary$fraction_cases,
  mis2_lower = apply(fraction_cases_mis2_boot, 1, quantile, 0.025),
  mis2_upper = apply(fraction_cases_mis2_boot, 1, quantile, 0.975),
  combined = ptv_summary$fraction_cases * autism_data$ptv_scale_factor +
    mis2_summary$fraction_cases,
  combined_lower = apply(fraction_cases_boot, 1, quantile, 0.025),
  combined_upper = apply(fraction_cases_boot, 1, quantile, 0.975)
)

variant_classes <- c("PTV", "Mis2", "Mis1", "Mis0", "Syn")
variant_model_names <- paste("Combined Probands", variant_classes)
mutational_variance <- do.call(rbind, lapply(seq_along(variant_classes), function(i) {
  model <- get_model(variant_model_names[i])
  scaling <- if (variant_classes[i] == "PTV") autism_data$ptv_scale_factor else 1
  data.frame(
    variant_class = variant_classes[i],
    estimate = model$mutvar_output$total_mutvar * scaling,
    lower = model$mutvar_output$mutvar_CI[1] * scaling,
    upper = model$mutvar_output$mutvar_CI[2] * scaling
  )
}))

ptv_enrichment <- data.frame(
  annotation = c("LOEUF1_mu1", "LOEUF1_mu2", paste0("LOEUF", 2:5)),
  estimate = ptv_model$mutvar_output$enrichment,
  lower = ptv_model$mutvar_output$enrich_CI[1, ],
  upper = ptv_model$mutvar_output$enrich_CI[2, ],
  fraction_mutvar = ptv_model$mutvar_output$frac_mutvar,
  fraction_expected = ptv_model$mutvar_output$frac_expected
)

penetrance_classes <- c("PTV", "Mis2", "Mis1")
effective_penetrance <- do.call(rbind, lapply(penetrance_classes, function(class) {
  model <- get_model(paste("Combined Probands", class))
  data.frame(
    variant_class = class,
    estimate = model$penetrance$effective_penetrance,
    lower = model$penetrance$effective_penetrance_CI[1],
    upper = model$penetrance$effective_penetrance_CI[2]
  )
}))

gene_mutational_variance <- data.frame(
  gene_id = names(combined_mutvar_gene),
  mutational_variance = as.numeric(combined_mutvar_gene),
  fraction_mutational_variance = combined_mutvar_gene / sum(combined_mutvar_gene),
  cumulative_fraction = cumsum(combined_mutvar_gene) / sum(combined_mutvar_gene)
)

summary_output <- list(
  metadata = list(
    model_file = normalizePath(model_file),
    optimizer = ptv_model$optimizer,
    prevalence = prevalence,
    prevalence_index = prevalence_index,
    ptv_scale_factor = autism_data$ptv_scale_factor,
    n_boot = n_boot,
    generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  ),
  mutational_variance = mutational_variance,
  cohort_table = cohort_table,
  polygenicity = polygenicity,
  half_mutvar_gene_count = half_index,
  enrichment = ptv_enrichment,
  fraction_cases = fraction_cases,
  effective_penetrance = effective_penetrance,
  gene_mutational_variance = gene_mutational_variance,
  ptv_number_genes = data.frame(
    threshold = thresholds,
    estimate = ptv_summary$number_genes,
    lower = apply(number_genes_boot, 1, quantile, 0.025),
    upper = apply(number_genes_boot, 1, quantile, 0.975)
  ),
  ptv_gene_rr_probability = ptv_gene_rr_probability
)

saveRDS(summary_output, file.path(output_dir, "figure2_summary.rds"))
write.table(
  mutational_variance,
  file.path(output_dir, "mutational_variance.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  cohort_table,
  file.path(output_dir, "figure2A_cohort_table.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  fraction_cases,
  file.path(output_dir, "fraction_cases_by_rr.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  effective_penetrance,
  file.path(output_dir, "effective_penetrance.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  gene_mutational_variance,
  file.path(output_dir, "gene_mutational_variance.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("Figure 2B-F summaries written to", output_dir, "\n")
cat("Genes explaining half of PTV+Mis2 mutational variance:", half_index, "\n")
print(mutational_variance, row.names = FALSE)
