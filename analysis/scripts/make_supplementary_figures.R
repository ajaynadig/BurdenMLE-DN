# One public entry point for all numbered supplementary figures that are not
# emitted alongside a main figure. Individual files remain small, readable
# implementation modules; study reproduction should call this script.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_file <- normalizePath(script_file, mustWork = TRUE)
example_dir <- dirname(script_file)
module_dir <- file.path(example_dir, "figure_modules")
repo_dir <- dirname(dirname(example_dir))
final_runs_dir <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT", unset = file.path(repo_dir, "analysis")), mustWork = FALSE)

modules <- c(
  "make_secondary_autism_figures.R",
  "make_sensitivity_figures.R",
  "make_posterior_diagnostic_figures.R",
  "make_ddid_figures.R",
  "make_autism_dd_prevalence_figure.R"
)
for (module in modules) {
  message("Creating supplementary figures via ", module)
  sys.source(file.path(module_dir, module), envir = new.env(parent = globalenv()))
}

# The no-overlap sensitivity deliberately reuses the Figure 4 plotting code.
rscript <- file.path(R.home("bin"), "Rscript")
status <- system2(
  rscript,
  c(
    file.path(example_dir, "make_figure4_autism_dd.R"),
    "--summary-file", file.path(
      final_runs_dir, "outputs", "derived", "sensitivities", "no_overlap_summary.rds"
    ),
    "--figure-dir", file.path(final_runs_dir, "outputs", "figures", "supplementary"),
    "--output-file", file.path(
      final_runs_dir, "outputs", "figures", "supplementary",
      "SupplementaryFigure15_AutismDDNoOverlap.pdf"
    ),
    "--panel-title", shQuote("Overlapping Samples Removed from Autism Dataset"),
    "--width", "14"
  )
)
if (status != 0) stop("No-overlap supplementary figure failed with status ", status)

cat("Numbered supplementary figures written to outputs/figures/supplementary.\n")
