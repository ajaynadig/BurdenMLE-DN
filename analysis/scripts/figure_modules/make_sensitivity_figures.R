library(ggplot2)
library(patchwork)
library(SummarizedExperiment)

# Sensitivity-figure entry point for the CES-gene exclusion analysis.

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
data_dir <- file.path(final_runs_dir, "outputs", "data")
figure_dir <- get_arg(
  "--figure-dir", file.path(final_runs_dir, "outputs", "figures", "supplementary")
)
derived_dir <- get_arg(
  "--derived-dir", file.path(final_runs_dir, "outputs", "derived", "sensitivities")
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)

for (source_file in c(
  "BurdenMLE_DN.R", "EM.R", "estimate_mutvar.R", "io.R",
  "likelihoods.R", "model.R", "secondary_analysis_functions.R"
)) source(if (identical(source_file, "secondary_analysis_functions.R")) file.path(repo_dir, "analysis", "scripts", source_file) else file.path(repo_dir, "R", source_file))

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
palette_variantclass <- c(
  "PTV" = "#FF5F8A", "Mis2" = "darkorange", "Mis1" = "#F9B332",
  "Mis0" = "#EFD09F", "Syn" = "grey50"
)
latest_model <- function(pattern, explicit, label) {
  candidates <- sort(list.files(data_dir, pattern = pattern, full.names = TRUE))
  path <- get_arg(explicit, if (length(candidates)) tail(candidates, 1) else NA_character_)
  if (is.na(path) || !file.exists(path)) stop("No ", label, " model output found.")
  path
}

# CES-gene exclusion sensitivity: original two-panel plot.
no_ces_file <- latest_model(
  "^models_autism_noCES_.*\\.Rdata$", "--no-ces-model-file", "no-CES"
)
no_ces_env <- new.env(parent = globalenv())
load(no_ces_file, envir = no_ces_env)
no_ces_table <- mutvar_enrichment_table(
  no_ces_env$autism_data, no_ces_env$BurdenMLE_DN_models_autism_noCES
)$mutvar
no_ces_plot_data <- no_ces_table[
  no_ces_table$prev_factor == no_ces_env$autism_data$prev_factors[1] &
    grepl("Combined", no_ces_table$name),
]
make_no_ces_panel <- function(role) {
  ggplot(
    no_ces_plot_data[no_ces_plot_data$role == role, ],
    aes(x = factor(variant_class, levels = c("Syn", "Mis0", "Mis1", "Mis2", "PTV")),
        y = mutvar, ymin = mutvar_lower, ymax = mutvar_upper, fill = variant_class)
  ) +
    geom_hline(yintercept = 0) +
    geom_pointrange(shape = 21, color = "black", size = 1) +
    scale_fill_manual(values = palette_variantclass) +
    theme_bhr_legend_gridlines() +
    labs(x = "Variant class", y = "Mutational variance\nObserved scale") +
    guides(fill = "none") +
    coord_cartesian(ylim = c(0, max(no_ces_plot_data$mutvar_upper, na.rm = TRUE) * 1.05))
}
no_ces_figure <-
  (make_no_ces_panel("Proband") + ggtitle("Probands")) +
  (make_no_ces_panel("Sibling") + ggtitle("Unaffected siblings")) +
  plot_annotation(tag_levels = "A")
ggsave(file.path(figure_dir, "SupplementaryFigure3_CESGeneExclusion.pdf"), no_ces_figure,
       device = cairo_pdf, width = 8, height = 3.7)
write.table(no_ces_plot_data, file.path(derived_dir, "no_ces_mutational_variance.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Sensitivity figures written to", figure_dir, "\n")
