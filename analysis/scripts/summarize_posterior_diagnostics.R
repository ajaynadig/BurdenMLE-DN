library(SummarizedExperiment)

# TADA comparison and coding-sequence-length summaries for the current PTV fit.

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
reference_file <- get_arg(
  "--reference-file",
  file.path(final_runs_dir, "inputs", "reference", "full_results_ASD_all_NPDs_2026-03-04.txt")
)
constraint_file <- get_arg(
  "--constraint-file",
  file.path(final_runs_dir, "inputs", "constraint", "gnomad.v4.1.constraint_metrics.tsv")
)
output_dir <- get_arg(
  "--output-dir", file.path(final_runs_dir, "outputs", "derived", "main_autism")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
for (path in c(model_file, reference_file, constraint_file)) {
  if (is.na(path) || !file.exists(path)) stop("Required input not found: ", path)
}
for (source_file in c("BurdenMLE_DN.R", "io.R", "model.R",
                      "secondary_analysis_functions.R")) {
  source(if (identical(source_file, "secondary_analysis_functions.R")) file.path(repo_dir, "analysis", "scripts", source_file) else file.path(repo_dir, "R", source_file))
}
load(model_file)
model_index <- setNames(seq_along(autism_data$loop_vars$names), autism_data$loop_vars$names)
ptv_index <- model_index[["Combined Probands PTV"]]
if (is.na(ptv_index)) stop("Combined Probands PTV model is missing.")
ptv_model <- BurdenMLE_DN_models_autism[[1]][[ptv_index]]
ptv_data <- get_genetic_data(ptv_index, autism_data)$genetic_data

# Match the external TADA results explicitly by Ensembl gene ID.
tada <- read.delim(reference_file, check.names = FALSE)
posterior <- as.data.frame(ptv_model$posterior_gene_estimates)
posterior$Gene_ID <- rownames(ptv_model$posterior_gene_estimates)
tada_keep <- tada[, c("Gene", "Gene_ID", "BF_PTV", "BF_Overall", "FDR", "Flag",
                      "p_PTV", "PTV_Proband")]
tada_comparison <- merge(posterior, tada_keep, by = "Gene_ID", all = FALSE, sort = FALSE)
tada_comparison <- tada_comparison[is.finite(tada_comparison$BF_PTV) &
                                     tada_comparison$BF_PTV > 0, ]
tada_stats <- data.frame(
  statistic = c("Spearman: posterior mean vs BF_PTV",
                "Spearman: -log10 posterior P(RR<=1) vs BF_PTV"),
  estimate = c(
    cor(tada_comparison$Posterior_Mean, tada_comparison$BF_PTV,
        method = "spearman", use = "complete.obs"),
    cor(-log10(tada_comparison$Posterior_ProbLessEqualZero), tada_comparison$BF_PTV,
        method = "spearman", use = "complete.obs")
  )
)

# gnomAD v4 contains gene-level and transcript-level rows. Select the canonical
# Ensembl transcript rows before matching; matching the first row for each gene
# instead selects the Entrez gene-level row, where CDS length is missing.
constraint <- read.delim(constraint_file, check.names = FALSE)
constraint <- constraint[
  is.finite(constraint$cds_length) & constraint$canonical == "true" &
    grepl("^ENSG", constraint$gene_id),
]
constraint <- constraint[!duplicated(constraint$gene_id), ]
gene_ids <- rowData(autism_data$counts)$gene_id
cds_length <- constraint$cds_length[match(gene_ids, constraint$gene_id)]
if (mean(is.finite(cds_length)) < 0.9) {
  stop("Fewer than 90% of fitted genes matched a gnomAD v4 CDS length.")
}
cds_length[!is.finite(cds_length)] <- median(cds_length, na.rm = TRUE)
thresholds <- 0:10

posterior_probability_grid <- function(model, genetic_data) {
  component_posterior <- (model$features %*% model$delta) * model$conditional_likelihood
  component_posterior <- component_posterior / rowSums(component_posterior)
  mu_grid <- seq(0.05, 1, by = 1 / model$grid_size)
  probability <- matrix(0, nrow(genetic_data),
                        length(model$component_endpoints) * length(mu_grid))
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
mean_cds_by_threshold <- function(model, genetic_data, cds) {
  grid <- posterior_probability_grid(model, genetic_data)
  indicators <- vapply(thresholds, function(x) as.numeric(exp(grid$log_rr) > x),
                       numeric(length(grid$log_rr)))
  gene_probabilities <- grid$probability %*% indicators
  as.numeric(crossprod(cds, gene_probabilities) / colSums(gene_probabilities))
}

mean_cds <- mean_cds_by_threshold(ptv_model, ptv_data, cds_length)
n_boot <- length(ptv_model$bootstrap_output$bootstrap_delta)
mean_cds_boot <- matrix(NA_real_, length(thresholds), n_boot)
progress <- txtProgressBar(min = 0, max = n_boot, style = 3)
for (iteration in seq_len(n_boot)) {
  sample_indices <- ptv_model$bootstrap_output$bootstrap_samples[, iteration]
  boot_model <- ptv_model
  boot_model$conditional_likelihood <- ptv_model$conditional_likelihood[sample_indices, , drop = FALSE]
  boot_model$features <- ptv_model$features[sample_indices, , drop = FALSE]
  boot_model$delta <- ptv_model$bootstrap_output$bootstrap_delta[[iteration]]
  mean_cds_boot[, iteration] <- mean_cds_by_threshold(
    boot_model, ptv_data[sample_indices, , drop = FALSE], cds_length[sample_indices]
  )
  setTxtProgressBar(progress, iteration)
}
close(progress)
cds_summary <- data.frame(
  threshold = thresholds, estimate = mean_cds,
  lower = apply(mean_cds_boot, 1, quantile, 0.025, na.rm = TRUE),
  upper = apply(mean_cds_boot, 1, quantile, 0.975, na.rm = TRUE)
)

write.table(tada_comparison, file.path(output_dir, "tada_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tada_stats, file.path(output_dir, "tada_comparison_statistics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cds_summary, file.path(output_dir, "cds_by_rr.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(tada = tada_comparison, tada_statistics = tada_stats, cds = cds_summary),
        file.path(output_dir, "posterior_diagnostics_summary.rds"))
cat("Posterior diagnostics written to", output_dir, "\n")
print(tada_stats, row.names = FALSE)
print(cds_summary, row.names = FALSE)
