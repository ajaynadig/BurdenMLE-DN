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

stopifnot(all(is.finite(results$tada$BF_PTV)), all(results$tada$BF_PTV > 0),
          all(is.finite(results$tada$Posterior_Mean)), all(results$tada$Posterior_Mean > 0))
xr <- range(log10(results$tada$BF_PTV))
yr <- range(log10(results$tada$Posterior_Mean))
xb <- seq(xr[1] - .001, xr[2] + .001, length.out = 61)
yb <- seq(yr[1] - .001, yr[2] + .001, length.out = 41)
comparison_plot <- function(data, title) {
  rho <- cor(data$Posterior_Mean, data$BF_PTV, method = "spearman")
  ggplot(data, aes(x = BF_PTV, y = Posterior_Mean)) +
    stat_bin_2d(breaks = list(x = 10^xb, y = 10^yb)) +
    scale_fill_gradient(low = "#d9e5ef", high = "#153d62", trans = "log10",
                        limits = c(1, nrow(results$tada)), breaks = c(1,10,100,1000,10000),
                        name = "Genes per bin",
                        guide = guide_colourbar(barwidth = grid::unit(7, "cm"))) +
    scale_y_log10(limits = 10^range(yb), breaks = c(1,3,10,30)) +
    scale_x_log10(limits = 10^range(xb)) +
    labs(title = title,
         subtitle = sprintf("n = %s; Spearman ρ = %.3f",
                            format(nrow(data), big.mark = ","), rho),
         x = "TADA Bayes Factor\n(PTV)",
         y = "BurdenMLE-DN Posterior Mean Rate Ratio\n(PTV)") +
    theme_bw(base_size = 12) +
    theme(axis.line = element_line(colour = "black"),
          axis.text = element_text(color = "black"),
          plot.title = element_text(size = 13), plot.tag = element_text(size = 22))
}
nominal_ptv <- !is.na(results$tada$p_PTV) & results$tada$p_PTV < 0.05
a <- comparison_plot(results$tada, "All genes: binned density")
b <- comparison_plot(results$tada[nominal_ptv, ], "Nominal PTV P < 0.05: binned density")
stopifnot(sum(ggplot_build(a)$data[[1]]$count) == nrow(results$tada),
          sum(ggplot_build(b)$data[[1]]$count) == sum(nominal_ptv))
BF_PosteriorMean_compare <- (a + b) + plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")
ggsave(file.path(figure_dir, "SupplementaryFigure6_TADAComparison.pdf"), BF_PosteriorMean_compare,
       device = cairo_pdf, width = 11, height = 5.5, dpi = "retina")
ggsave(file.path(figure_dir, "SupplementaryFigure6_TADAComparison.png"), BF_PosteriorMean_compare,
       width = 11, height = 5.5, dpi = 180)

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
                                format(length_summary$n, big.mark = ",", trim = TRUE))
GeneSetCDS_plot <- ggplot(
  length_summary,
  aes(x = factor(label, levels = label), y = median)
) +
  geom_col(width = 0.5, color = "black", fill = "white") +
  theme_bhr_gridlines() +
  theme(axis.text.x = element_text(size = 11)) +
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
