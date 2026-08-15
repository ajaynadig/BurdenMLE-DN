cli_candidates <- c(
  file.path("analysis", "scripts", "run_models_cli.R"),
  file.path("..", "analysis", "scripts", "run_models_cli.R")
)
cli_file <- cli_candidates[file.exists(cli_candidates)][1L]

# analysis/ is intentionally excluded from package builds, so this focused test
# runs in the source checkout and is a no-op in an installed-package check.
if (length(cli_file) && !is.na(cli_file)) {
  repo_dir <- normalizePath(
    file.path(dirname(cli_file), "..", ".."), mustWork = TRUE
  )
  source(cli_file)
  defaults <- parse_run_models_args(character(), as.Date("2026-08-15"))
  explicit_em <- parse_run_models_args(c("--optimizer", "EM"))
  explicit_mixsqp <- parse_run_models_args("--optimizer=MixSQP")
  stopifnot(
    identical(defaults$optimizer, "mixsqp"),
    identical(defaults$run_date, "Aug15_26"),
    identical(explicit_em$optimizer, "EM"),
    identical(explicit_mixsqp$optimizer, "mixsqp"),
    grepl("default: mixsqp", run_models_help(), fixed = TRUE),
    isTRUE(parse_run_models_args("--help")$help)
  )

  expect_error <- function(args) {
    result <- try(parse_run_models_args(args), silent = TRUE)
    stopifnot(inherits(result, "try-error"))
  }
  expect_error(c("--optimizer", "other"))
  expect_error(c("--optimizer"))

  help_output <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(file.path(repo_dir, "analysis", "scripts", "run_models.R"), "--help"),
    stdout = TRUE,
    stderr = TRUE
  )
  stopifnot(any(grepl("default: mixsqp", help_output, fixed = TRUE)))

  reproduction <- readLines(
    file.path(repo_dir, "analysis", "reproduce_study.sh"), warn = FALSE
  )
  stopifnot(
    any(grepl("--optimizer EM", reproduction, fixed = TRUE)),
    any(grepl("--optimizer mixsqp", reproduction, fixed = TRUE))
  )
}
