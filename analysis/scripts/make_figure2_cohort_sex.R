library(ggplot2)
library(patchwork)

# Add the original cohort-by-sex panels G-H to the manuscript-style Figure 2.

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
example_dir <- dirname(script_file)
repo_dir <- dirname(example_dir)
final_runs_dir <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT", unset = file.path(repo_dir, "analysis")), mustWork = FALSE)
summary_file <- get_arg(
  "--cohort-summary-file",
  file.path(final_runs_dir, "outputs", "derived", "cohort_ddid", "cohort_sex_summary.rds")
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
if (!file.exists(summary_file)) stop("Cohort/sex summary not found: ", summary_file)
dataset_mutvar_results <- readRDS(summary_file)

# Creates cohort_table_panel and panels B-F using the manuscript ggplot code.
source(file.path(example_dir, "make_figure2_main_autism.R"))

dataset_mutvar_maindisplay <- ggplot(
  data = dataset_mutvar_results[dataset_mutvar_results$sex != "Both", ],
  mapping = aes(
    x = factor(study, levels = c("SPARK", "ASC", "GeneDX")),
    y = mutvar_combined, ymin = mutvar_combined_lower, ymax = mutvar_combined_upper,
    fill = sex
  )
) +
  geom_hline(yintercept = 0) +
  geom_pointrange(
    position = position_dodge(width = 0.5), color = "black", shape = 21, size = 1
  ) +
  theme_bhr_legend_gridlines() +
  scale_fill_manual(values = c("Male" = "chartreuse3", "Female" = "#9467bd")) +
  labs(x = "Study", y = "Mutational Variance\nPTV + Mis2", color = "Sex") +
  guides(fill = "none")

dataset_fraccase_maindisplay <- ggplot(
  data = dataset_mutvar_results[dataset_mutvar_results$sex != "Both", ],
  mapping = aes(
    x = factor(study, levels = c("SPARK", "ASC", "GeneDX")),
    y = fraccase_RR5, ymin = fraccase_RR5_lower, ymax = fraccase_RR5_upper,
    fill = sex
  )
) +
  geom_hline(yintercept = 0) +
  geom_col(color = "black", position = position_dodge(width = 1), width = 0.8) +
  geom_errorbar(position = position_dodge(width = 1), width = 0.1) +
  theme_bhr_legend_gridlines() +
  scale_fill_manual(values = c("Male" = "chartreuse3", "Female" = "#9467bd")) +
  labs(x = "Study", y = "Fraction of Cases\nRR>5, PTV + Mis2", fill = "Sex")

combined_figure2 <-
  (cohort_table_panel | Proband_Vmu_plot) /
  (Polygenicity_plot | EnrichmentPlot) /
  (FracCases_plot | peneff_plot) /
  (dataset_fraccase_maindisplay | dataset_mutvar_maindisplay) +
  plot_annotation(tag_levels = "A", tag_prefix = "", tag_suffix = "") &
  theme(plot.tag = element_text(size = 25))

ggsave(file.path(figure_dir, "Figure2.pdf"), combined_figure2,
       device = cairo_pdf, width = 12, height = 15)

# Cohort-only panels retained as a supplementary/reviewer-facing output.
dataset_mutvar_revision <- ggplot(
  data = dataset_mutvar_results[dataset_mutvar_results$sex == "Both", ],
  mapping = aes(
    x = factor(study, levels = c("SPARK", "ASC", "GeneDX")),
    y = mutvar_combined, ymin = mutvar_combined_lower, ymax = mutvar_combined_upper
  )
) +
  geom_hline(yintercept = 0) +
  geom_pointrange(position = position_dodge(width = 0.5), color = "black", size = 1) +
  theme_bhr_legend_gridlines() +
  labs(x = "Study", y = "Mutational Variance\nPTV + Mis2", color = "Sex") +
  guides(fill = "none")

dataset_fraccase_revision <- ggplot(
  data = dataset_mutvar_results[dataset_mutvar_results$sex == "Both", ],
  mapping = aes(
    x = factor(study, levels = c("SPARK", "ASC", "GeneDX")),
    y = fraccase_RR5, ymin = fraccase_RR5_lower, ymax = fraccase_RR5_upper
  )
) +
  geom_hline(yintercept = 0) +
  geom_col(color = "black", position = position_dodge(width = 1),
           width = 0.8, fill = "white") +
  geom_errorbar(position = position_dodge(width = 1), width = 0.1) +
  theme_bhr_legend_gridlines() +
  labs(x = "Study", y = "Fraction of Cases\nRR>5, PTV + Mis2", fill = "Sex")

cohort_only_figure <- dataset_mutvar_revision | dataset_fraccase_revision
ggsave(file.path(supp_figure_dir, "SupplementaryFigure10_CohortHeterogeneity.pdf"),
       cohort_only_figure, device = cairo_pdf, width = 8, height = 4)
write.table(dataset_mutvar_results, file.path(table_dir, "cohort_sex_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Full manuscript-style Figure 2 written to", figure_dir, "\n")
