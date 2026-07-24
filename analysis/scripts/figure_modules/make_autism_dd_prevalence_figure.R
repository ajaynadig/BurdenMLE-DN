library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_file <- normalizePath(script_file, mustWork = TRUE)
repo_dir <- dirname(dirname(script_file))
final_runs_dir <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT", unset = file.path(repo_dir, "analysis")), mustWork = FALSE)
input_file <- file.path(final_runs_dir, "outputs", "derived", "autism_dd", "prevalence_sweep.tsv")
figure_dir <- file.path(final_runs_dir, "outputs", "figures", "supplementary")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
summary <- read.delim(input_file)
theme_bhr_legend_gridlines <- function() {
  theme_bw() + theme(axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 15, color = "black"),
    axis.title = element_text(size = 15, color = "black"),
    legend.text = element_text(size = 15), legend.title = element_text(size = 15))
}
plot <- ggplot(summary, aes(x = prevalence, y = estimate, ymin = lower, ymax = upper,
                            color = dataset, fill = dataset)) +
  scale_color_manual(values = c("Autism" = "darkgreen", "DD" = "orchid2")) +
  scale_fill_manual(values = c("Autism" = "forestgreen", "DD" = "orchid1")) +
  geom_hline(yintercept = 0) + geom_line() +
  geom_ribbon(alpha = 0.2, color = "grey40") +
  labs(x = "Prevalence", y = "Mutational Variance\nPTV+Mis2", fill = "Diagnosis") +
  theme_bhr_legend_gridlines() + guides(color = "none")
ggsave(file.path(figure_dir, "SupplementaryFigure14_DDPrevalenceSensitivity.pdf"), plot,
       device = cairo_pdf, width = 6, height = 4, dpi = "retina")
cat("Autism/DD prevalence figure written to", figure_dir, "\n")
