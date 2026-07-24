library(ggplot2)
library(patchwork)
library(grid)

# This script intentionally preserves the original manuscript ggplot design.
# Refactoring here is limited to paths, data inputs, and dynamic annotations.

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
summary_file <- get_arg(
  "--summary-file",
  file.path(final_runs_dir, "outputs", "derived", "main_autism", "figure2_summary.rds")
)
figure_dir <- get_arg(
  "--figure-dir",
  file.path(final_runs_dir, "outputs", "figures", "main")
)
table_dir <- get_arg(
  "--table-dir",
  file.path(final_runs_dir, "outputs", "tables", "manuscript_values")
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(summary_file)) stop("Summary file does not exist: ", summary_file)
results <- readRDS(summary_file)

theme_bhr_legend_gridlines <- function() {
  theme_bw() +
    theme(
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size = 15, color = "black"),
      axis.title = element_text(size = 15, color = "black"),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 15),
      strip.text.x = element_text(size = 15),
      strip.background = element_rect(fill = "white")
    )
}

palette_variantclass <- c(
  "PTV" = "#FF5F8A", "Mis2" = "darkorange", "Mis1" = "#F9B332",
  "Mis0" = "#EFD09F", "Syn" = "grey50"
)

# Figure 2A: draw the manuscript table as grid primitives so the values remain
# automated while retaining the original sparse-rule design.
make_cohort_table_grob <- function(cohort_table) {
  display <- cohort_table
  display$n_proband <- format(display$n_proband, big.mark = ",", scientific = FALSE)
  display$n_constrained_ptv <- format(
    display$n_constrained_ptv, big.mark = ",", scientific = FALSE
  )
  display$rr_constrained_ptv <- sprintf("%.1f", display$rr_constrained_ptv)

  x <- c(0.11, 0.40, 0.67, 0.91)
  y_header <- 0.88
  y_rows <- seq(0.70, 0.08, length.out = nrow(display))
  all_separator_y <- mean(tail(y_rows, 2))
  grobs <- list(
    textGrob(expression(italic(Dataset)), x = x[1], y = y_header,
             gp = gpar(fontsize = 12.5)),
    textGrob(expression(italic(N)[Proband]), x = x[2], y = y_header,
             gp = gpar(fontsize = 12.5)),
    textGrob(expression(italic(N)[Cons-PTV]), x = x[3], y = y_header,
             gp = gpar(fontsize = 12.5)),
    textGrob(expression(italic(RR)[Cons-PTV]), x = x[4], y = y_header,
             gp = gpar(fontsize = 12.5)),
    linesGrob(x = c(0.02, 0.98), y = c(0.79, 0.79), gp = gpar(lwd = 1)),
    linesGrob(x = c(0.02, 0.98), y = rep(all_separator_y, 2), gp = gpar(lwd = 1)),
    linesGrob(x = c(0.27, 0.27), y = c(0.02, 0.79), gp = gpar(lwd = 0.8)),
    linesGrob(x = c(0.55, 0.55), y = c(0.02, 0.79), gp = gpar(lwd = 0.8)),
    linesGrob(x = c(0.80, 0.80), y = c(0.02, 0.79), gp = gpar(lwd = 0.8))
  )
  for (i in seq_len(nrow(display))) {
    fontface <- if (display$dataset[i] == "All") "bold" else "plain"
    values <- c(
      display$dataset[i], display$n_proband[i],
      display$n_constrained_ptv[i], display$rr_constrained_ptv[i]
    )
    for (j in seq_along(values)) {
      grobs[[length(grobs) + 1]] <- textGrob(
        values[j], x = x[j], y = y_rows[i],
        gp = gpar(fontsize = 14, fontface = fontface)
      )
    }
  }
  do.call(grobTree, grobs)
}

cohort_table_grob <- make_cohort_table_grob(results$cohort_table)
cohort_table_panel <- wrap_elements(full = cohort_table_grob)

# Figure 2B: retain the original vertical pointrange design.
plot_data <- results$mutational_variance
Proband_Vmu_plot <- ggplot(
  data = plot_data,
  mapping = aes(
    x = factor(variant_class, levels = c("Syn", "Mis0", "Mis1", "Mis2", "PTV")),
    y = estimate, ymin = lower, ymax = upper,
    fill = factor(variant_class, levels = c("Syn", "Mis0", "Mis1", "Mis2", "PTV"))
  )
) +
  geom_hline(yintercept = 0) +
  geom_pointrange(position = position_dodge2(width = 0.25), shape = 21,
                  color = "black", size = 1) +
  scale_fill_manual(values = palette_variantclass) +
  theme_bhr_legend_gridlines() +
  labs(x = "Variant Class", y = "Mutational Variance\nObserved Scale") +
  theme(
    strip.text.y = element_text(size = 15, angle = 270),
    axis.text.x = element_text(color = rev(palette_variantclass), face = "bold")
  ) +
  ylim(0, 0.045) +
  guides(fill = "none")

# Figure 2C: retain absolute cumulative mutational variance on the y-axis.
poly <- results$polygenicity
Polygenicity_plot <- ggplot(
  data = poly,
  mapping = aes(x = number_genes, y = cumulative_mutvar_absolute)
) +
  geom_ribbon(
    aes(ymin = lower_absolute, ymax = upper_absolute),
    alpha = 0.2, color = "grey40", fill = "deepskyblue3"
  ) +
  geom_line(color = "black") +
  theme_bhr_legend_gridlines() +
  geom_vline(xintercept = results$half_mutvar_gene_count, linetype = "dashed") +
  scale_x_log10() +
  labs(x = "Number of Genes", y = "Cumulative Mutational Variance\nPTV + Mis2") +
  annotate(
    "text", x = results$half_mutvar_gene_count + 10, y = 0.005,
    label = paste0(results$half_mutvar_gene_count,
                   " genes explain half\nof mutational variance"),
    size = 5, fontface = "italic", hjust = 0
  )

# Figure 2D: retain the original bars, enrichment error bars, and fraction labels.
enrichment_viz_df <- results$enrichment
enrichment_viz_df$annot_reformat <- gsub("_", "\n", enrichment_viz_df$annotation)
enrichment_viz_df$variant_class <- "PTV"
EnrichmentPlot <- ggplot(
  enrichment_viz_df[enrichment_viz_df$annotation != "LOEUF5", ],
  mapping = aes(
    x = factor(annot_reformat, levels = annot_reformat[annotation != "LOEUF5"]),
    y = estimate, ymin = lower, ymax = upper,
    fill = factor(variant_class, levels = "PTV"),
    label = paste0(round(fraction_mutvar * 100), "%")
  )
) +
  geom_col(color = "black") +
  geom_errorbar(color = "black", width = 0.2) +
  geom_label(
    aes(y = -0.25), fill = "white", fontface = "bold", size = 4.5,
    hjust = 0.4, label.r = unit(0, "lines")
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = palette_variantclass) +
  theme_bhr_legend_gridlines() +
  labs(x = "Annotation", y = "Enrichment\nMutational Variance", fill = "Variant Class") +
  guides(fill = "none") +
  theme(axis.text.x = element_text(size = 15))

# Figure 2E: retain the original line/ribbon and ACMG annotation.
FracCases_plot <- ggplot(
  data = results$fraction_cases,
  mapping = aes(
    x = threshold, y = combined,
    ymin = combined_lower, ymax = combined_upper
  )
) +
  geom_hline(yintercept = 0) +
  geom_ribbon(alpha = 0.2, color = "grey40", fill = "deepskyblue3") +
  geom_line() +
  labs(x = "RR Threshold", y = "Fraction of Cases\nHighly Penetrant PTV/Mis2") +
  theme_bhr_legend_gridlines() +
  geom_vline(xintercept = 5, linetype = "dashed") +
  ylim(0, 0.105) +
  annotate(
    "text", x = 5.5, y = 0.095, label = "ACMG\nCriterion",
    size = 5, fontface = "italic", hjust = 0
  )

# Figure 2F: retain the original vertical effective-penetrance design.
peneff_df <- results$effective_penetrance
peneff_plot <- ggplot(
  data = peneff_df,
  mapping = aes(
    x = factor(variant_class, levels = c("Mis1", "Mis2", "PTV")),
    y = estimate, ymin = lower, ymax = upper,
    fill = factor(variant_class, levels = c("Mis1", "Mis2", "PTV"))
  )
) +
  geom_hline(yintercept = 0) +
  geom_pointrange(shape = 21, color = "black", size = 1) +
  labs(x = "Variant Class", y = "Effective Penetrance") +
  scale_fill_manual(values = palette_variantclass) +
  theme_bhr_legend_gridlines() +
  theme(axis.text.x = element_text(color = rev(palette_variantclass)[3:5], face = "bold")) +
  ylim(0, 1) +
  guides(fill = "none")

# These panel objects are assembled with G-H by make_figure2_cohort_sex.R.
# Individual main-text subpanel PDFs are intentionally not written.

manuscript_values <- rbind(
  data.frame(
    quantity = paste("mutational_variance", results$mutational_variance$variant_class),
    estimate = results$mutational_variance$estimate,
    lower = results$mutational_variance$lower,
    upper = results$mutational_variance$upper
  ),
  data.frame(
    quantity = paste("effective_penetrance", results$effective_penetrance$variant_class),
    estimate = results$effective_penetrance$estimate,
    lower = results$effective_penetrance$lower,
    upper = results$effective_penetrance$upper
  )
)
write.table(
  manuscript_values, file.path(table_dir, "figure2_main_values.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
cat("Prepared manuscript-style Figure 2A-F panel objects for final assembly.\n")
