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
repo_dir <- get0("repo_dir", envir = parent.env(environment()), inherits = TRUE)
if (is.null(repo_dir)) {
  script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  script_file <- normalizePath(script_file, mustWork = TRUE)
  repo_dir <- dirname(dirname(dirname(dirname(script_file))))
}
final_runs_dir <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT", unset = file.path(repo_dir, "analysis")), mustWork = FALSE)
summary_file <- get_arg(
  "--summary-file",
  file.path(final_runs_dir, "outputs", "derived", "main_autism", "posterior_diagnostics_summary.rds")
)
figure_dir <- get_arg(
  "--figure-dir", file.path(final_runs_dir, "outputs", "figures", "supplementary")
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(summary_file)) stop("Posterior diagnostics summary not found.")
results <- readRDS(summary_file)
theme_bhr_gridlines <- function() {
  theme_bw() + theme(
    axis.line = element_line(colour = "black"), axis.text = element_text(size = 15, color = "black"),
    axis.title = element_text(size = 15, color = "black"), legend.position = "none"
  )
}

BF_PosteriorMean_compare <- ggplot(
  results$tada, aes(x = BF_PTV, y = Posterior_Mean)
) +
  geom_point() + scale_y_log10() + scale_x_log10() +
  labs(x = "TADA Bayes Factor\n(PTV)",
       y = "BurdenMLE-DN Posterior Mean Rate Ratio\n(PTV)") +
  theme_bhr_gridlines()
ggsave(file.path(figure_dir, "SupplementaryFigure6_TADAComparison.pdf"), BF_PosteriorMean_compare,
       device = cairo_pdf, width = 5, height = 5, dpi = "retina")

# Coding-length comparison for the GER and neuronal-communication gene sets.
# This was Supplementary Figure 5 in the original submission and remains a
# separate analysis from the posterior-mean CDS plot below.
constraint <- read.delim(file.path(
  final_runs_dir, "inputs", "constraint", "gnomad.v2.1.1.lof_metrics.by_gene.txt"
))
ger <- read.delim(file.path(final_runs_dir, "inputs", "reference", "GER_IDs.txt"),
                  header = FALSE)[[1]]
nc <- read.delim(file.path(final_runs_dir, "inputs", "reference", "NC_IDs.txt"),
                 header = FALSE)[[1]]
length_table <- data.frame(
  gene = constraint$gene_id,
  length = constraint$cds_length,
  group = "Neither"
)
length_table$group[length_table$gene %in% ger] <- "GER"
length_table$group[length_table$gene %in% nc] <- "NC"
groups <- c("Neither", "GER", "NC")
length_summary <- data.frame(
  group = groups,
  median = vapply(groups, function(group) {
    median(length_table$length[length_table$group == group], na.rm = TRUE)
  }, numeric(1)),
  n = vapply(groups, function(group) sum(length_table$group == group), integer(1))
)
length_summary$label <- sprintf("%s\n(N = %s)", length_summary$group,
                                format(length_summary$n, big.mark = ","))
GeneSetCDS_plot <- ggplot(
  length_summary,
  aes(x = factor(label, levels = label), y = median)
) +
  geom_col(width = 0.5, color = "black", fill = "white") +
  theme_bhr_gridlines() +
  labs(x = "Gene Set", y = "Coding Length\nMedian")
ggsave(file.path(figure_dir, "SupplementaryFigure8_GeneSetCodingLength.pdf"),
       GeneSetCDS_plot, device = cairo_pdf, width = 5, height = 4, dpi = "retina")

MeanCDS_plot <- ggplot(results$cds, aes(x = threshold, y = estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2,
              color = "grey40", fill = "deepskyblue3") +
  geom_line(color = "black") + theme_bhr_gridlines() +
  labs(x = "Rate Ratio Threshold", y = "Mean Coding Sequence Length")
ggsave(file.path(figure_dir, "SupplementaryFigure9_CodingLengthByRateRatio.pdf"), MeanCDS_plot,
       device = cairo_pdf, width = 8, height = 4, dpi = "retina")
cat("Posterior diagnostic figures written to", figure_dir, "\n")
