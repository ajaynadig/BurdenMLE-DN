library(ggplot2)
library(patchwork)

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
  file.path(final_runs_dir, "outputs", "derived", "cohort_ddid", "ddid_summary.rds")
)
figure_dir <- get_arg(
  "--figure-dir", file.path(final_runs_dir, "outputs", "figures", "supplementary")
)
table_dir <- get_arg(
  "--table-dir", file.path(final_runs_dir, "outputs", "tables", "supplementary")
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(summary_file)) stop("DDID summary not found: ", summary_file)
DDID_compare_df <- readRDS(summary_file)

theme_bhr_legend_gridlines <- function() {
  theme_bw() + theme(
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 15, color = "black"),
    axis.title = element_text(size = 15, color = "black"),
    legend.text = element_text(size = 15), legend.title = element_text(size = 15),
    strip.text.x = element_text(size = 15), strip.text.y = element_text(size = 15),
    strip.background = element_rect(fill = "white")
  )
}

combined <- DDID_compare_df[DDID_compare_df$Dataset == "Combined", ]
DDID_mutvar_plot <- ggplot(
  data = combined,
  mapping = aes(x = DDID, y = mutvar_combined,
                ymin = mutvar_combined_lower, ymax = mutvar_combined_upper)
) +
  geom_hline(mapping = aes(yintercept = 0)) + geom_pointrange() +
  theme_bhr_legend_gridlines() + ylim(0, 0.045) +
  labs(x = "", y = "Mutational Variance\nPTV + Mis2")

DDID_fraccase_plot <- ggplot(
  data = combined,
  mapping = aes(x = DDID, y = fraccase_RR5,
                ymin = fraccase_RR5_lower, ymax = fraccase_RR5_upper)
) +
  geom_hline(mapping = aes(yintercept = 0)) + geom_pointrange() +
  theme_bhr_legend_gridlines() +
  labs(x = "", y = "Fraction of Cases\nRR>5, PTV + Mis2")

DDID_effRR_plot <- ggplot(
  data = combined,
  mapping = aes(x = DDID, y = PTV_effRR,
                ymin = PTV_effRR_lower, ymax = PTV_effRR_upper)
) +
  geom_hline(mapping = aes(yintercept = 0)) + geom_pointrange() +
  theme_bhr_legend_gridlines() +
  labs(x = "", y = "Effective Rate Ratio\nPTV")

DDID_fig <- DDID_mutvar_plot + DDID_fraccase_plot + DDID_effRR_plot +
  plot_annotation(tag_levels = "A", tag_prefix = "", tag_suffix = "") &
  theme(plot.tag = element_text(size = 25))

maxstrat_table <- DDID_compare_df[DDID_compare_df$Dataset != "Combined", ]
DDID_mutvar_plot_maxstrat <- ggplot(
  data = maxstrat_table,
  mapping = aes(x = DDID, y = mutvar_combined,
                ymin = mutvar_combined_lower, ymax = mutvar_combined_upper)
) +
  geom_hline(mapping = aes(yintercept = 0)) + geom_pointrange() +
  theme_bhr_legend_gridlines() +
  labs(x = "", y = "Mutational Variance\nPTV + Mis2") +
  facet_grid(rows = vars(Sex), cols = vars(Dataset))

DDID_fraccase_plot_maxstrat <- ggplot(
  data = maxstrat_table,
  mapping = aes(x = DDID, y = fraccase_RR5,
                ymin = fraccase_RR5_lower, ymax = fraccase_RR5_upper)
) +
  geom_hline(mapping = aes(yintercept = 0)) + geom_pointrange() +
  theme_bhr_legend_gridlines() +
  labs(x = "", y = "Fraction of Cases\nRR>5, PTV + Mis2") +
  facet_grid(rows = vars(Sex), cols = vars(Dataset))

DDID_effRR_plot_maxstrat <- ggplot(
  data = maxstrat_table,
  mapping = aes(x = DDID, y = PTV_effRR,
                ymin = PTV_effRR_lower, ymax = PTV_effRR_upper)
) +
  geom_hline(mapping = aes(yintercept = 0)) + geom_pointrange() +
  theme_bhr_legend_gridlines() +
  labs(x = "", y = "Effective Rate Ratio\nPTV") +
  facet_grid(rows = vars(Sex), cols = vars(Dataset))

maxstrat_figure <- DDID_mutvar_plot_maxstrat /
  DDID_fraccase_plot_maxstrat /
  DDID_effRR_plot_maxstrat +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 25))

ggsave(file.path(figure_dir, "SupplementaryFigure11_DDIDComparison.pdf"),
       DDID_fig, device = cairo_pdf, width = 12, height = 4)
ggsave(file.path(figure_dir, "SupplementaryFigure12_DDIDMaximallyStratified.pdf"),
       maxstrat_figure, device = cairo_pdf, width = 12, height = 20)
write.table(DDID_compare_df, file.path(table_dir, "ddid_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("DDID supplementary figures written to", figure_dir, "\n")
