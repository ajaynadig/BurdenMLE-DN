repo_candidates <- c(".", "..")
repo_dir <- repo_candidates[dir.exists(file.path(repo_candidates, "analysis", "scripts"))][1L]

# analysis/ is intentionally excluded from package builds. This test runs in
# the source checkout and is a no-op inside the built package check.
if (length(repo_dir) && !is.na(repo_dir)) {
  repo_dir <- normalizePath(repo_dir, mustWork = TRUE)
  scripts_dir <- file.path(repo_dir, "analysis", "scripts")
  reproduction <- readLines(
    file.path(repo_dir, "analysis", "reproduce_study.sh"), warn = FALSE
  )

  matches <- regmatches(
    reproduction,
    gregexpr('\\$\\{SCRIPT_DIR\\}/[A-Za-z0-9_.-]+\\.R', reproduction)
  )
  stages <- unique(sub('^\\$\\{SCRIPT_DIR\\}/', '', unlist(matches)))
  stopifnot(
    length(stages) > 0L,
    all(file.exists(file.path(scripts_dir, stages)))
  )

  r_files <- list.files(
    scripts_dir, pattern = "\\.R$", recursive = TRUE, full.names = TRUE
  )
  invisible(lapply(r_files, parse))

  top_level <- list.files(scripts_dir, pattern = "\\.R$", full.names = TRUE)
  top_text <- unlist(lapply(top_level, readLines, warn = FALSE), use.names = FALSE)
  stopifnot(
    !any(grepl("repo_dir <- dirname\\(dirname\\(script_file\\)\\)", top_text)),
    !any(grepl("repo_dir <- dirname\\(example_dir\\)", top_text))
  )

  dependencies <- c(
    file.path(repo_dir, "R", c(
      "BurdenMLE_DN.R", "EM.R", "MixSQP.R", "estimate_mutvar.R",
      "io.R", "likelihoods.R", "model.R"
    )),
    file.path(scripts_dir, c(
      "secondary_analysis_functions.R", "model_artifacts.R", "run_models_cli.R",
      "set_up_asc.R", "set_up_asc_NoKaplanisOverlap.R", "set_up_kaplanis.R"
    ))
  )
  stopifnot(all(file.exists(dependencies)))
}

