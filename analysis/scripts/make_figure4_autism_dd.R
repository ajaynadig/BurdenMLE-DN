library(ggplot2)
library(patchwork)

# Manuscript Figure 4: autism versus developmental-disorder estimates.
# The plotting code deliberately retains the design from visualization_trio.R.

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
summary_file <- get_arg(
  "--summary-file",
  file.path(final_runs_dir, "outputs", "derived", "autism_dd", "figure4_summary.rds")
)
figure_dir <- get_arg(
  "--figure-dir", file.path(final_runs_dir, "outputs", "figures", "main")
)
output_file <- get_arg("--output-file", file.path(figure_dir, "Figure4.pdf"))
panel_title <- get_arg("--panel-title", "")
figure_width <- as.numeric(get_arg("--width", "11"))
figure_height <- as.numeric(get_arg("--height", "11"))
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(summary_file)) stop("Figure 4 summary not found: ", summary_file)
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

mutvar_compare_plot <- ggplot(
  data = results$mutational_variance,
  mapping = aes(
    x = factor(variant_class, levels = c("Syn", "Mis0", "Mis1", "Mis2", "PTV")),
    y = estimate, ymin = lower, ymax = upper, fill = diagnosis
  )
) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = c("Autism" = "darkgreen", "DD" = "orchid2")) +
  geom_pointrange(position = position_dodge2(width = 0.4), shape = 21, size = 1) +
  theme_bhr_legend_gridlines() +
  labs(x = "Variant Class", y = "Mutational Variance\nObserved Scale") +
  theme(
    strip.text.y = element_text(size = 15, angle = 270),
    axis.text.x = element_text(color = rev(palette_variantclass), face = "bold")
  ) +
  ylim(0, 0.1) +
  guides(fill = "none") +
  geom_rect(mapping = aes(xmin = 0.85, xmax = 2.4, ymin = 0.030, ymax = 0.044),
            color = "black", fill = "white", inherit.aes = FALSE) +
  annotate("text", x = 1, y = 0.04, label = "Autism", color = "darkgreen",
           size = 6, fontface = "bold", hjust = 0) +
  annotate("text", x = 1, y = 0.034, label = "DD", color = "orchid3",
           size = 6, fontface = "bold", hjust = 0)

FracCases_compare_plot <- ggplot(
  data = results$fraction_cases,
  mapping = aes(
    x = threshold, y = estimate, ymin = lower, ymax = upper,
    color = diagnosis, fill = diagnosis
  )
) +
  scale_color_manual(values = c("Autism" = "darkgreen", "DD" = "orchid2")) +
  scale_fill_manual(values = c("Autism" = "forestgreen", "DD" = "orchid1")) +
  geom_hline(yintercept = 0) +
  geom_line() +
  geom_ribbon(alpha = 0.2, color = "grey40") +
  labs(x = "RR Threshold", y = "Fraction of Cases\nHighly Penetrant PTV/Mis2") +
  theme_bhr_legend_gridlines() +
  ylim(0, 0.22) +
  geom_vline(xintercept = 5, linetype = "dashed") +
  guides(color = "none", fill = "none") +
  annotate("text", x = 5.5, y = 0.21, label = "ACMG\nCriterion",
           size = 5, fontface = "italic", hjust = 0)

RReff_compare_plot <- ggplot(
  data = results$effective_rr,
  mapping = aes(
    x = factor(variant_class, levels = c("Mis1", "Mis2", "PTV")),
    y = estimate, ymin = lower, ymax = upper, fill = diagnosis
  )
) +
  scale_fill_manual(values = c("Autism" = "darkgreen", "DD" = "orchid2")) +
  geom_hline(yintercept = 0) +
  geom_pointrange(position = position_dodge2(width = 0.4), shape = 21, size = 1) +
  labs(x = "Variant Class", y = "Effective Rate Ratio") +
  theme_bhr_legend_gridlines() +
  theme(axis.text.x = element_text(color = rev(palette_variantclass)[-c(1, 2)],
                                   face = "bold")) +
  # The current DD PTV upper confidence limit is slightly above 80; using an
  # 85-unit limit keeps the complete interval instead of dropping its segment.
  ylim(0, 85) +
  guides(fill = "none")

numgenes_compare_plot <- ggplot(
  data = results$number_genes,
  mapping = aes(
    x = threshold, y = estimate, ymin = lower, ymax = upper,
    color = diagnosis, fill = diagnosis
  )
) +
  scale_color_manual(values = c("Autism" = "darkgreen", "DD" = "orchid2")) +
  scale_fill_manual(values = c("Autism" = "forestgreen", "DD" = "orchid1")) +
  geom_hline(yintercept = 0) +
  geom_line() +
  geom_ribbon(alpha = 0.2, color = "grey40") +
  labs(x = "RR Threshold", y = "Number of Genes with\nPTV Effect Size > Threshold") +
  theme_bhr_legend_gridlines() +
  guides(color = "none", fill = "none") +
  scale_y_log10()

if (nzchar(panel_title)) {
  mutvar_compare_plot <- mutvar_compare_plot + ggtitle(panel_title)
  FracCases_compare_plot <- FracCases_compare_plot + ggtitle(panel_title)
  RReff_compare_plot <- RReff_compare_plot + ggtitle(panel_title)
  numgenes_compare_plot <- numgenes_compare_plot + ggtitle(panel_title)
}

figure4 <- mutvar_compare_plot + FracCases_compare_plot +
  RReff_compare_plot + numgenes_compare_plot +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20, face = "bold"))

ggsave(output_file, figure4,
       device = cairo_pdf, width = figure_width, height = figure_height, dpi = "retina")
cat("Full manuscript-style autism/DD figure written to", output_file, "\n")
