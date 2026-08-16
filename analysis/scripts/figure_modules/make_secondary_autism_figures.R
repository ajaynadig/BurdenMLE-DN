library(ggplot2)

# Preserve the supplementary ggplot designs from visualization_trio.R.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit)) return(sub(paste0("^", flag, "="), "", hit[1]))
  i <- match(flag, args)
  if (!is.na(i)) return(args[i + 1])
  default
}
repo_dir <- get0("repo_dir", envir = parent.env(environment()), inherits = TRUE)
if (is.null(repo_dir)) {
  script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  script_file <- normalizePath(script_file, mustWork = TRUE)
  repo_dir <- dirname(dirname(dirname(dirname(script_file))))
}
final_runs_dir <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT", unset = file.path(repo_dir, "analysis")), mustWork = FALSE)
summary_file <- get_arg(
  "--summary-file",
  file.path(final_runs_dir, "outputs", "derived", "main_autism", "secondary_autism_summary.rds")
)
figure_dir <- get_arg(
  "--figure-dir", file.path(final_runs_dir, "outputs", "figures", "supplementary")
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(summary_file)) stop("Secondary autism summary not found.")
results <- readRDS(summary_file)

theme_bhr_legend_gridlines <- function() {
  theme_bw() + theme(
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 15, color = "black"),
    axis.title = element_text(size = 15, color = "black"),
    legend.text = element_text(size = 15), legend.title = element_text(size = 15),
    strip.text.x = element_text(size = 15), strip.background = element_rect(fill = "white")
  )
}
palette_variantclass <- c(
  "PTV" = "#FF5F8A", "Mis2" = "darkorange", "Mis1" = "#F9B332",
  "Mis0" = "#EFD09F", "Syn" = "grey50"
)

sibling_data <- results$sibling_mutational_variance
Sibling_Vmu_plot <- ggplot(
  sibling_data,
  aes(x = factor(variant_class, levels = c("Syn", "Mis0", "Mis1", "Mis2", "PTV")),
      y = estimate, ymin = lower, ymax = upper,
      fill = factor(variant_class, levels = c("Syn", "Mis0", "Mis1", "Mis2", "PTV")))
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
  ylim(0, 0.0375) + guides(fill = "none") +
  ggtitle(sprintf("Unaffected Siblings; N = %s",
                  format(unique(sibling_data$n), big.mark = ",")))
ggsave(file.path(figure_dir, "SupplementaryFigure2_SiblingNegativeControl.pdf"), Sibling_Vmu_plot,
       device = cairo_pdf, width = 4.5, height = 3, dpi = "retina")

prevalence_data <- results$prevalence_mutational_variance
prevalence_data$prevalence_label <- paste0(
  format(round(prevalence_data$prevalence * 100, 2), trim = TRUE), "%"
)
prevalence_levels <- prevalence_data$prevalence_label[
  order(prevalence_data$prevalence)
]
prevalence_levels <- unique(prevalence_levels)
Prevalence_MutVar_plot <- ggplot(
  prevalence_data,
  aes(x = factor(variant_class, levels = c("Syn", "Mis0", "Mis1", "Mis2", "PTV")),
      y = estimate, ymin = lower, ymax = upper,
      fill = factor(prevalence_label, levels = prevalence_levels))
) +
  geom_hline(yintercept = 0) +
  geom_pointrange(position = position_dodge2(width = 0.6), shape = 21,
                  color = "black", size = 1) +
  theme_bhr_legend_gridlines() +
  labs(x = "Variant Class", y = "Mutational Variance\nObserved Scale", fill = "Prevalence") +
  theme(
    strip.text.y = element_text(size = 15, angle = 270),
    axis.text.x = element_text(color = rev(palette_variantclass), face = "bold")
  )
ggsave(file.path(figure_dir, "SupplementaryFigure5_PrevalenceSensitivity.pdf"), Prevalence_MutVar_plot,
       device = cairo_pdf, width = 6, height = 4, dpi = "retina")
cat("Sibling and prevalence figures written to", figure_dir, "\n")
