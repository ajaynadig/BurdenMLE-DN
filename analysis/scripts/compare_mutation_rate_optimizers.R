# Compare EM and MixSQP for the gnomAD synonymous mutation-rate calibration.
# This is a diagnostic only and does not alter any manuscript model objects.

library(ggplot2)
library(patchwork)

script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(script_file)) {
  script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
}
script_file <- normalizePath(script_file, mustWork = TRUE)
repo_dir <- dirname(dirname(script_file))
final_runs_dir <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT", unset = file.path(repo_dir, "analysis")), mustWork = FALSE)

for (source_file in c(
  "BurdenMLE_DN.R", "EM.R", "MixSQP.R", "estimate_mutvar.R", "io.R",
  "likelihoods.R", "model.R", "secondary_analysis_functions.R"
)) {
  source(if (identical(source_file, "secondary_analysis_functions.R")) file.path(repo_dir, "analysis", "scripts", source_file) else file.path(repo_dir, "R", source_file))
}

gnomad <- read.table(
  file.path(final_runs_dir, "inputs", "constraint", "gnomad.v4.1.constraint_metrics.tsv"),
  header = TRUE
)
keep <- gnomad$canonical == "true" & grepl("ENSG", gnomad$gene_id)
calibration_data <- data.frame(
  case_count = gnomad$syn.obs[keep],
  expected_count = gnomad$syn.exp[keep],
  gene_id = gnomad$gene_id[keep]
)
calibration_data <- calibration_data[complete.cases(calibration_data), ]
rownames(calibration_data) <- calibration_data$gene_id

fit_calibration <- function(optimizer, max_iter = 1000, tol = 1e-6) {
  elapsed <- system.time({
    fit <- BurdenMLE_DN(
      calibration_data,
      features = NULL,
      component_endpoints = seq(-2, 2, length.out = 31),
      mutvar_est = FALSE,
      max_iter = max_iter,
      tol = tol,
      null_sim = FALSE,
      bootstrap = FALSE,
      return_likelihood = TRUE,
      estimate_posteriors = TRUE,
      estimate_effective_penetrance = FALSE,
      optimizer = optimizer
    )
  })[["elapsed"]]
  list(fit = fit, elapsed = elapsed)
}

em <- fit_calibration("EM")
mixsqp <- fit_calibration("mixsqp")
# This tighter EM fit diagnoses whether any difference from MixSQP is caused by
# the production EM stopping rule rather than a different likelihood target.
em_tight <- fit_calibration("EM", max_iter = 10000, tol = 1e-10)

em_posteriors <- transform(
  em$fit$posterior_gene_estimates,
  gene_id = calibration_data$gene_id,
  em_factor = Posterior_Mean
)
mixsqp_posteriors <- transform(
  mixsqp$fit$posterior_gene_estimates,
  gene_id = calibration_data$gene_id,
  mixsqp_factor = Posterior_Mean
)
tight_em_factor <- em_tight$fit$posterior_gene_estimates$Posterior_Mean
posterior_comparison <- merge(
  em_posteriors[, c("gene_id", "em_factor")],
  mixsqp_posteriors[, c("gene_id", "mixsqp_factor")],
  by = "gene_id"
)
posterior_comparison$difference <-
  posterior_comparison$mixsqp_factor - posterior_comparison$em_factor
posterior_comparison$tight_em_factor <- tight_em_factor
posterior_comparison$mixsqp_tight_em_difference <-
  posterior_comparison$mixsqp_factor - posterior_comparison$tight_em_factor

weight_comparison <- data.frame(
  component = em$fit$component_endpoints,
  em_weight = drop(em$fit$delta),
  mixsqp_weight = drop(mixsqp$fit$delta)
)

em_ll <- tail(em$fit$ll[is.finite(em$fit$ll)], 1)
mixsqp_ll <- as.numeric(mixsqp$fit$ll)
tight_em_ll <- tail(em_tight$fit$ll[is.finite(em_tight$fit$ll)], 1)
diagnostics <- data.frame(
  metric = c(
    "genes", "EM log likelihood", "MixSQP log likelihood",
    "MixSQP - EM log likelihood", "maximum absolute factor difference",
    "median absolute factor difference", "factor correlation",
    "EM iterations", "EM elapsed seconds", "MixSQP elapsed seconds",
    "Tight EM log likelihood", "MixSQP - tight EM log likelihood",
    "tight EM iterations",
    "maximum absolute MixSQP - tight EM factor difference",
    "median absolute MixSQP - tight EM factor difference",
    "MixSQP - tight EM factor correlation"
  ),
  value = c(
    nrow(posterior_comparison), em_ll, mixsqp_ll, mixsqp_ll - em_ll,
    max(abs(posterior_comparison$difference)),
    median(abs(posterior_comparison$difference)),
    cor(posterior_comparison$em_factor, posterior_comparison$mixsqp_factor),
    sum(is.finite(em$fit$ll)), em$elapsed, mixsqp$elapsed,
    tight_em_ll, mixsqp_ll - tight_em_ll, sum(is.finite(em_tight$fit$ll)),
    max(abs(posterior_comparison$mixsqp_tight_em_difference)),
    median(abs(posterior_comparison$mixsqp_tight_em_difference)),
    cor(posterior_comparison$mixsqp_factor, posterior_comparison$tight_em_factor)
  )
)

factor_plot <- ggplot(
  posterior_comparison, aes(em_factor, mixsqp_factor)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.25, size = 0.7, color = "#345995") +
  coord_equal() +
  theme_bw() +
  labs(
    x = "EM posterior mutation-rate correction factor",
    y = "MixSQP posterior mutation-rate correction factor"
  )

weight_plot <- ggplot(
  weight_comparison, aes(em_weight, mixsqp_weight)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(size = 2, color = "#E45756") +
  coord_equal() +
  theme_bw() +
  labs(x = "EM mixture weight", y = "MixSQP mixture weight")

diagnostic_figure <- factor_plot + weight_plot +
  plot_annotation(title = "Mutation-rate calibration: MixSQP versus EM")

figure_dir <- file.path(final_runs_dir, "outputs", "figures", "diagnostics", "optimizer")
table_dir <- file.path(final_runs_dir, "outputs", "tables", "diagnostics", "optimizer")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  file.path(figure_dir, "MutationRateCalibrationOptimizerComparison.pdf"),
  diagnostic_figure, device = cairo_pdf, width = 10, height = 4.8
)
write.table(
  posterior_comparison,
  file.path(table_dir, "mutation_rate_calibration_gene_factors.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  weight_comparison,
  file.path(table_dir, "mutation_rate_calibration_mixture_weights.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  diagnostics,
  file.path(table_dir, "mutation_rate_calibration_optimizer_diagnostics.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

print(diagnostics, row.names = FALSE)
