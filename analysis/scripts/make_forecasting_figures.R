library(ggplot2)
library(patchwork)
library(MetBrewer)

# Preserve the revised manuscript forecasting ggplot code. This script only
# replaces hardcoded paths and model indices with derived, name-based inputs.

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

main_summary_file <- get_arg(
  "--main-summary-file",
  file.path(final_runs_dir, "outputs", "derived", "main_autism", "figure2_summary.rds")
)
forecast_file <- get_arg(
  "--forecast-file",
  file.path(final_runs_dir, "outputs", "derived", "forecasting", "forecast_summary.rds")
)
number_genes_file <- get_arg(
  "--number-genes-file",
  file.path(final_runs_dir, "outputs", "derived", "autism_dd", "number_genes.tsv")
)
reference_file <- get_arg(
  "--reference-file",
  file.path(final_runs_dir, "inputs", "reference", "full_results_ASD_all_NPDs_2026-03-04.txt")
)
figure_dir <- get_arg(
  "--figure-dir", file.path(final_runs_dir, "outputs", "figures", "main")
)
supp_figure_dir <- get_arg(
  "--supp-figure-dir", file.path(final_runs_dir, "outputs", "figures", "supplementary")
)
table_dir <- get_arg(
  "--table-dir", file.path(final_runs_dir, "outputs", "tables", "supplementary")
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

for (file in c(main_summary_file, forecast_file, number_genes_file, reference_file)) {
  if (!file.exists(file)) stop("Required input does not exist: ", file)
}
main_summary <- readRDS(main_summary_file)
forecast <- readRDS(forecast_file)
number_genes <- read.table(
  number_genes_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE
)
significance <- read.table(reference_file, sep = "\t", header = TRUE,
                           stringsAsFactors = FALSE, check.names = FALSE)

theme_bhr_legend_gridlines <- function() {
  theme_bw() + theme(
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 15, color = "black"),
    axis.title = element_text(size = 15, color = "black"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    strip.text.x = element_text(size = 15),
    strip.background = element_rect(fill = "white")
  )
}
theme_bhr_gridlines <- function() {
  theme_bhr_legend_gridlines() + theme(legend.position = "none")
}

thresholds <- c(2, 5, 10, 20)
gene_probability <- main_summary$ptv_gene_rr_probability
probability_columns <- paste0("rr_greater_", thresholds)
if (!all(probability_columns %in% names(gene_probability))) {
  stop("Main summary predates Figure 3 gene probabilities; rerun summarize_main_autism.R.")
}

prediction <- data.frame(
  RR_thresh = thresholds,
  num_genes = colSums(gene_probability[probability_columns])
)
prediction_intervals <- number_genes[
  number_genes$diagnosis == "Autism" & number_genes$threshold %in% thresholds,
  c("threshold", "estimate", "lower", "upper")
]
prediction_intervals <- prediction_intervals[
  match(thresholds, prediction_intervals$threshold), , drop = FALSE
]
if (anyNA(prediction_intervals) ||
    !isTRUE(all.equal(prediction$num_genes, prediction_intervals$estimate,
                      tolerance = 1e-10))) {
  stop("Figure 3A bootstrap intervals do not match its point estimates.")
}
prediction$lower <- prediction_intervals$lower
prediction$upper <- prediction_intervals$upper
significant_gene_counts <- function(fdr_threshold) {
  genes <- significance$Gene_ID[
    !significance$Flag & !is.na(significance$FDR) & significance$FDR < fdr_threshold
  ]
  matched <- gene_probability$gene_id %in% genes
  data.frame(
    RR_thresh = thresholds,
    num_genes = colSums(gene_probability[matched, probability_columns, drop = FALSE]),
    lower = NA_real_,
    upper = NA_real_
  )
}
asc_fdr_005 <- significant_gene_counts(0.05)
asc_fdr_0001 <- significant_gene_counts(0.001)

# Revised response-letter visualization: grouped bars, not a display table.
plot_df <- rbind(
  data.frame(prediction, study = "BurdenMLE-DN\nPrediction"),
  data.frame(asc_fdr_005, study = "ASC 2026\nFDR < 0.05"),
  data.frame(asc_fdr_0001, study = "ASC 2026\nFDR < 0.001")
)
plot_df$RR_thresh <- factor(
  plot_df$RR_thresh,
  levels = c(2, 5, 10, 20),
  labels = c("\u2265 2", "\u2265 5", "\u2265 10", "\u2265 20")
)
plot_df$study <- factor(
  plot_df$study,
  levels = c("BurdenMLE-DN\nPrediction", "ASC 2026\nFDR < 0.05", "ASC 2026\nFDR < 0.001")
)
plot_df$bar_label <- format(
  round(plot_df$num_genes), big.mark = ",", scientific = FALSE, trim = TRUE
)
plot_df$label_y <- ifelse(is.na(plot_df$upper), plot_df$num_genes, plot_df$upper)
dodge_width <- 0.82

p_effectsize_thresholds <- ggplot(
  plot_df,
  aes(x = RR_thresh, y = num_genes, fill = study, group = study)
) +
  geom_col(
    position = position_dodge(width = dodge_width),
    width = 0.72, colour = "black"
  ) +
  geom_errorbar(
    data = plot_df,
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = dodge_width),
    width = 0.08, linewidth = 0.7
  ) +
  geom_text(
    aes(y = label_y, label = bar_label),
    position = position_dodge(width = dodge_width),
    vjust = -0.45, size = 5, fontface = "bold"
  ) +
  scale_fill_manual(
    name = "Gene set",
    values = c(
      "BurdenMLE-DN\nPrediction" = "#3D334A",
      "ASC 2026\nFDR < 0.05" = "#9BC4B8",
      "ASC 2026\nFDR < 0.001" = "#3F7C70"
    ),
    breaks = c(
      "BurdenMLE-DN\nPrediction", "ASC 2026\nFDR < 0.05", "ASC 2026\nFDR < 0.001"
    )
  ) +
  scale_y_continuous(
    labels = function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE),
    expand = expansion(mult = c(0, 0.11))
  ) +
  labs(x = "Risk-ratio threshold", y = "Number of significant genes") +
  theme_bhr_gridlines() +
  theme(
    legend.position = "right", legend.direction = "vertical",
    legend.justification = "center", panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(margin = margin(t = 6)
  ))

forecast_df <- forecast$cohort
palette <- MetBrewer::met.brewer("Juarez", 3, "discrete")
ForecastingPlot <- ggplot(
  data = forecast_df[forecast_df$threshold == 0, ],
  mapping = aes(
    x = total_sample_size, y = count_bonferroni_mean,
    ymin = count_bonferroni_mean - count_bonferroni_sd,
    ymax = count_bonferroni_mean + count_bonferroni_sd,
    color = factor(scenario, levels = c("SPARK", "ASC", "GeneDx")),
    fill = factor(scenario, levels = c("SPARK", "ASC", "GeneDx"))
  )
) +
  geom_ribbon(alpha = 0.2, linetype = "blank") +
  geom_line() +
  geom_line(
    mapping = aes(y = count_bonferroni_mean - count_bonferroni_sd), alpha = 0.3
  ) +
  geom_line(
    mapping = aes(y = count_bonferroni_mean + count_bonferroni_sd), alpha = 0.3
  ) +
  scale_color_manual(values = palette) +
  scale_fill_manual(values = palette) +
  theme_bhr_legend_gridlines() +
  labs(
    x = "Sample Size", y = "Number of Significant Genes\n(PTV)",
    fill = "New Data\nSource"
  ) +
  guides(alpha = "none", color = "none",
         fill = guide_legend(override.aes = list(alpha = 1)))

ForecastFig_All <- p_effectsize_thresholds / ForecastingPlot +
  plot_annotation(tag_levels = "A", tag_prefix = "", tag_suffix = "") &
  theme(plot.tag = element_text(size = 25))

ddid_df <- forecast$ddid
palette_DDID <- MetBrewer::met.brewer("Egypt", 2, "discrete")
ForecastingPlot_DDID <- ggplot(
  data = ddid_df[ddid_df$threshold == 0, ],
  mapping = aes(
    x = total_sample_size, y = count_bonferroni_mean,
    ymin = count_bonferroni_mean - count_bonferroni_sd,
    ymax = count_bonferroni_mean + count_bonferroni_sd,
    color = factor(scenario, levels = c("DDID", "Non-DDID")),
    fill = factor(scenario, levels = c("DDID", "Non-DDID"))
  )
) +
  geom_ribbon(alpha = 0.2, linetype = "blank") +
  geom_line() +
  geom_line(mapping = aes(y = count_bonferroni_mean - count_bonferroni_sd), alpha = 0.3) +
  geom_line(mapping = aes(y = count_bonferroni_mean + count_bonferroni_sd), alpha = 0.3) +
  scale_color_manual(values = palette_DDID) +
  scale_fill_manual(values = palette_DDID) +
  theme_bhr_legend_gridlines() +
  labs(
    x = "Sample Size", y = "Number of Significant Genes\n(PTV)",
    fill = "New Data\nSource"
  ) +
  guides(alpha = "none", color = "none",
         fill = guide_legend(override.aes = list(alpha = 1)))

ggsave(file.path(figure_dir, "Figure3.pdf"),
       ForecastFig_All, device = cairo_pdf, width = 9, height = 9)
ggsave(file.path(supp_figure_dir, "SupplementaryFigure13_DDIDForecasting.pdf"),
       ForecastingPlot_DDID, device = cairo_pdf, width = 9, height = 4)

write.table(
  plot_df[, c("RR_thresh", "num_genes", "lower", "upper", "study", "bar_label")],
  file.path(table_dir, "figure3_gene_count_bars.tsv"),
            sep = "\t", quote = TRUE, row.names = FALSE)
write.table(forecast_df[forecast_df$threshold == 0, ],
            file.path(table_dir, "forecasting_gene_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Forecasting figures written to", figure_dir, "\n")
print(plot_df[, c("RR_thresh", "study", "num_genes")], row.names = FALSE)
