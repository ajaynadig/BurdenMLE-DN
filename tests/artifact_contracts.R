library(pbapply)

for (source_file in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) {
  source(source_file)
}
source("analysis/scripts/model_artifacts.R")

expect_error <- function(code, pattern = NULL) {
  result <- tryCatch(
    {
      force(code)
      NULL
    },
    error = identity
  )
  stopifnot(inherits(result, "error"))
  if (!is.null(pattern)) stopifnot(grepl(pattern, conditionMessage(result)))
  invisible(result)
}

expect_warning <- function(code, pattern) {
  warning_message <- NULL
  value <- withCallingHandlers(
    force(code),
    warning = function(warning) {
      warning_message <<- conditionMessage(warning)
      invokeRestart("muffleWarning")
    }
  )
  stopifnot(!is.null(warning_message), grepl(pattern, warning_message))
  value
}

fixture_dir <- tempfile("artifact-contract-")
dir.create(fixture_dir)

gene_ids <- paste0("gene_", 1:4)
raw_input <- data.frame(
  case_count = c(3, 0, 1, 2),
  case_rate = c(0.002, 0.003, 0.004, 0.005),
  N = rep(100, 4),
  row.names = gene_ids
)
expected_input <- data.frame(
  case_count = raw_input$case_count,
  expected_count = 2 * raw_input$N * raw_input$case_rate,
  row.names = gene_ids
)
fit_args <- list(
  component_endpoints = c(0, log(5)),
  mutvar_est = FALSE,
  bootstrap = FALSE,
  null_sim = FALSE,
  estimate_posteriors = TRUE,
  estimate_effective_penetrance = FALSE,
  optimizer = "mixsqp"
)
raw_fit <- do.call(BurdenMLE_DN, c(list(input_data = raw_input), fit_args))
expected_fit <- do.call(
  BurdenMLE_DN, c(list(input_data = expected_input), fit_args)
)
stopifnot(
  isTRUE(raw_fit$fit_status$converged),
  isTRUE(expected_fit$fit_status$converged),
  identical(rownames(raw_fit$features), gene_ids),
  identical(rownames(expected_fit$features), gene_ids)
)

BurdenMLE_DN_models_autism_test <- list(list(raw_fit, expected_fit))
autism_data <- raw_input
artifact_path <- file.path(fixture_dir, "models_autism_test_original.Rdata")
save(BurdenMLE_DN_models_autism_test, autism_data, file = artifact_path)

manifest_path <- file.path(fixture_dir, "model_manifest_test_fixture.rds")
invisible(write_model_manifest(
  manifest_path,
  run_id = "fixture",
  mode = "test",
  repo_dir = ".",
  settings = list(
    optimizer = "mixsqp", seed = 1L, grid_size = 10L,
    bootstrap_count = 0L, max_effect_size = list(autism = 5)
  ),
  artifacts = c(test = artifact_path)
))

# A newer-looking file has no effect: the manifest-selected path is exact.
decoy_path <- file.path(fixture_dir, "models_autism_test_zzzz.Rdata")
decoy <- "not a model"
save(decoy, file = decoy_path)
loaded <- new.env(parent = globalenv())
info <- load_model_artifact(manifest_path, "test", envir = loaded)
manifest <- read_model_manifest(manifest_path)
stopifnot(
  identical(info$path, normalizePath(artifact_path)),
  identical(manifest$artifacts$test$gene_ids, gene_ids),
  identical(manifest$artifacts$test$fit_records$fit_id, c("1/1", "1/2")),
  identical(manifest$artifacts$test$fit_records$bootstrap_replicates,
            c(0L, 0L)),
  identical(loaded$BurdenMLE_DN_models_autism_test,
            BurdenMLE_DN_models_autism_test),
  identical(loaded$autism_data, autism_data)
)

# Package posterior calculations agree before and after the analysis loader.
direct_sampler <- posterior_gene_sampler(raw_fit, raw_input)
loaded_fit <- loaded$BurdenMLE_DN_models_autism_test[[1]][[1]]
loaded_sampler <- posterior_gene_sampler(loaded_fit, loaded$autism_data)
set.seed(90210)
direct_draws <- posterior_gene_samples(direct_sampler, 4L)
set.seed(90210)
loaded_draws <- posterior_gene_samples(loaded_sampler, 4L)
stopifnot(identical(direct_draws, loaded_draws))

# Evaluate the actual forecasting helper without executing its workflow.
forecast_expressions <- parse("analysis/scripts/forecasting_script_revision.R")
forecast_definition <- which(vapply(forecast_expressions, function(expression) {
  is.call(expression) && identical(expression[[1]], quote(`<-`)) &&
    identical(expression[[2]], as.name("forecast_once"))
}, logical(1)))
stopifnot(length(forecast_definition) == 1L)
forecast_environment <- new.env(parent = globalenv())
eval(forecast_expressions[[forecast_definition]], envir = forecast_environment)
set.seed(444)
direct_forecast <- forecast_environment$forecast_once(
  raw_input, direct_sampler, n_new = 25, gamma_scaling_factor = 0.8
)
set.seed(444)
loaded_forecast <- forecast_environment$forecast_once(
  loaded$autism_data, loaded_sampler, n_new = 25,
  gamma_scaling_factor = 0.8
)
stopifnot(identical(direct_forecast, loaded_forecast))

expect_error(load_model_artifact(NULL, "test"), "exactly one")
expect_error(
  load_model_artifact(manifest_path, "test", legacy_path = artifact_path),
  "exactly one"
)
expect_error(load_model_artifact(artifact_path, "test"))

bad_schema <- manifest
bad_schema$schema_version <- 99L
bad_schema_path <- file.path(fixture_dir, "bad_schema.rds")
saveRDS(bad_schema, bad_schema_path)
expect_error(read_model_manifest(bad_schema_path), "Unsupported")

absent_entry <- manifest
absent_entry$mode <- "main"
absent_entry$artifacts <- list(autism = manifest$artifacts$test)
absent_entry_path <- file.path(fixture_dir, "absent_entry.rds")
saveRDS(absent_entry, absent_entry_path)
expect_error(load_model_artifact(absent_entry_path, "ddd"), "does not contain")

missing_file <- manifest
missing_file$artifacts$test$path <- file.path(fixture_dir, "missing.Rdata")
missing_file_path <- file.path(fixture_dir, "missing_file.rds")
saveRDS(missing_file, missing_file_path)
expect_error(load_model_artifact(missing_file_path, "test"), "does not exist")

noncanonical <- manifest
dir.create(file.path(fixture_dir, "subdirectory"))
noncanonical$artifacts$test$path <- file.path(
  fixture_dir, "subdirectory", "..", basename(artifact_path)
)
noncanonical_path <- file.path(fixture_dir, "noncanonical.rds")
saveRDS(noncanonical, noncanonical_path)
expect_error(load_model_artifact(noncanonical_path, "test"), "not the exact")

wrong_gene_record <- manifest
wrong_gene_record$artifacts$test$gene_ids <- rev(gene_ids)
wrong_gene_path <- file.path(fixture_dir, "wrong_gene_record.rds")
saveRDS(wrong_gene_record, wrong_gene_path)
expect_error(load_model_artifact(wrong_gene_path, "test"), "gene universe")

wrong_status_record <- manifest
wrong_status_record$artifacts$test$fit_records$status_code[1] <- "other"
wrong_status_path <- file.path(fixture_dir, "wrong_status_record.rds")
saveRDS(wrong_status_record, wrong_status_path)
expect_error(load_model_artifact(wrong_status_path, "test"), "fit-status")

wrong_optimizer <- manifest
wrong_optimizer$settings$optimizer <- "EM"
wrong_optimizer_path <- file.path(fixture_dir, "wrong_optimizer.rds")
saveRDS(wrong_optimizer, wrong_optimizer_path)
expect_error(load_model_artifact(wrong_optimizer_path, "test"),
             "optimizer does not match")

bad_repository <- manifest
bad_repository$repository$commit <- NULL
bad_repository_path <- file.path(fixture_dir, "bad_repository.rds")
saveRDS(bad_repository, bad_repository_path)
expect_error(read_model_manifest(bad_repository_path),
             "repository metadata")

missing_object_path <- file.path(fixture_dir, "missing_object.Rdata")
save(BurdenMLE_DN_models_autism_test, file = missing_object_path)
missing_object_manifest <- manifest
missing_object_manifest$artifacts$test$path <- normalizePath(missing_object_path)
missing_object_manifest_path <- file.path(fixture_dir, "missing_object.rds")
saveRDS(missing_object_manifest, missing_object_manifest_path)
expect_error(load_model_artifact(missing_object_manifest_path, "test"),
             "missing required objects")

reordered_models <- BurdenMLE_DN_models_autism_test
reordered_models[[1]][[2]]$features <-
  reordered_models[[1]][[2]]$features[rev(gene_ids), , drop = FALSE]
reordered_models[[1]][[2]]$conditional_likelihood <-
  reordered_models[[1]][[2]]$conditional_likelihood[
    rev(gene_ids), , drop = FALSE
  ]
reordered_models[[1]][[2]]$conditional_log_likelihood <-
  reordered_models[[1]][[2]]$conditional_log_likelihood[
    rev(gene_ids), , drop = FALSE
  ]
reordered_models[[1]][[2]]$likelihood_log_scale <-
  reordered_models[[1]][[2]]$likelihood_log_scale[rev(gene_ids)]
BurdenMLE_DN_models_autism_test <- reordered_models
reordered_path <- file.path(fixture_dir, "reordered.Rdata")
save(BurdenMLE_DN_models_autism_test, autism_data, file = reordered_path)
expect_error(inspect_model_artifact(reordered_path, "test"),
             "gene universe")

misaligned_models <- list(list(raw_fit))
misaligned_models[[1]][[1]]$conditional_likelihood <-
  misaligned_models[[1]][[1]]$conditional_likelihood[
    rev(gene_ids), , drop = FALSE
  ]
BurdenMLE_DN_models_autism_test <- misaligned_models
misaligned_path <- file.path(fixture_dir, "misaligned_likelihood.Rdata")
save(BurdenMLE_DN_models_autism_test, autism_data, file = misaligned_path)
expect_error(inspect_model_artifact(misaligned_path, "test"),
             "conditional_likelihood")

missing_status_models <- list(list(raw_fit))
missing_status_models[[1]][[1]]$fit_status <- NULL
BurdenMLE_DN_models_autism_test <- missing_status_models
missing_status_path <- file.path(fixture_dir, "missing_status.Rdata")
save(BurdenMLE_DN_models_autism_test, autism_data, file = missing_status_path)
expect_error(inspect_model_artifact(missing_status_path, "test"),
             "required modern fit fields")

missing_likelihood_models <- list(list(raw_fit))
missing_likelihood_models[[1]][[1]]$conditional_likelihood <- NULL
BurdenMLE_DN_models_autism_test <- missing_likelihood_models
missing_likelihood_path <- file.path(fixture_dir, "missing_likelihood.Rdata")
save(BurdenMLE_DN_models_autism_test, autism_data,
     file = missing_likelihood_path)
expect_error(inspect_model_artifact(missing_likelihood_path, "test"),
             "required modern fit fields")

malformed_delta_models <- list(list(raw_fit))
malformed_delta_models[[1]][[1]]$delta[1, 1] <-
  malformed_delta_models[[1]][[1]]$delta[1, 1] + 0.1
BurdenMLE_DN_models_autism_test <- malformed_delta_models
malformed_delta_path <- file.path(fixture_dir, "malformed_delta.Rdata")
save(BurdenMLE_DN_models_autism_test, autism_data,
     file = malformed_delta_path)
expect_error(inspect_model_artifact(malformed_delta_path, "test"),
             "invalid mixture weights")

stale_status_models <- list(list(raw_fit))
stale_status_models[[1]][[1]]$delta <-
  stale_status_models[[1]][[1]]$delta[, 2:1, drop = FALSE]
BurdenMLE_DN_models_autism_test <- stale_status_models
stale_status_path <- file.path(fixture_dir, "stale_status_weights.Rdata")
save(BurdenMLE_DN_models_autism_test, autism_data, file = stale_status_path)
expect_error(inspect_model_artifact(stale_status_path, "test"),
             "inconsistent with fit_status")

unusable_models <- list(list(raw_fit))
unusable_models[[1]][[1]]$fit_status$usable <- FALSE
unusable_models[[1]][[1]]$fit_status$code <- "invalid_weights"
BurdenMLE_DN_models_autism_test <- unusable_models
unusable_path <- file.path(fixture_dir, "unusable.Rdata")
save(BurdenMLE_DN_models_autism_test, autism_data, file = unusable_path)
expect_error(inspect_model_artifact(unusable_path, "test"), "unusable")

nonconverged_models <- list(list(raw_fit))
nonconverged_models[[1]][[1]]$fit_status$converged <- FALSE
nonconverged_models[[1]][[1]]$fit_status$code <- "nonconverged"
BurdenMLE_DN_models_autism_test <- nonconverged_models
nonconverged_path <- file.path(fixture_dir, "nonconverged.Rdata")
save(BurdenMLE_DN_models_autism_test, autism_data, file = nonconverged_path)
expect_error(inspect_model_artifact(nonconverged_path, "test"),
             "did not converge")

unusable_bootstrap_models <- list(list(raw_fit))
unusable_bootstrap_status <- raw_fit$fit_status
unusable_bootstrap_status$usable <- FALSE
unusable_bootstrap_status$code <- "invalid_weights"
unusable_bootstrap_models[[1]][[1]]$bootstrap_output <- list(
  fit_status = list(bootstrap_1 = unusable_bootstrap_status),
  reliable = FALSE
)
unusable_bootstrap_models[[1]][[1]]$uncertainty_reliable <- FALSE
BurdenMLE_DN_models_autism_test <- unusable_bootstrap_models
unusable_bootstrap_path <- file.path(fixture_dir, "unusable_bootstrap.Rdata")
save(BurdenMLE_DN_models_autism_test, autism_data,
     file = unusable_bootstrap_path)
expect_error(inspect_model_artifact(unusable_bootstrap_path, "test"),
             "bootstrap/bootstrap_1.*unusable")

unreliable_models <- list(list(raw_fit))
nonconverged_bootstrap_status <- raw_fit$fit_status
nonconverged_bootstrap_status$converged <- FALSE
nonconverged_bootstrap_status$code <- "nonconverged"
unreliable_models[[1]][[1]]$bootstrap_output <- list(
  fit_status = list(bootstrap_1 = nonconverged_bootstrap_status),
  reliable = FALSE
)
unreliable_models[[1]][[1]]$uncertainty_reliable <- FALSE
BurdenMLE_DN_models_autism_test <- unreliable_models
unreliable_path <- file.path(fixture_dir, "unreliable.Rdata")
save(BurdenMLE_DN_models_autism_test, autism_data, file = unreliable_path)
unreliable_manifest <- file.path(fixture_dir, "unreliable_manifest.rds")
invisible(write_model_manifest(
  unreliable_manifest, "unreliable", "test", ".",
  settings = list(
    optimizer = "mixsqp", seed = 1L, grid_size = 10L,
    bootstrap_count = 1L, max_effect_size = list(autism = 5)
  ),
  artifacts = c(test = unreliable_path)
))
unreliable_info <- expect_warning(
  load_model_artifact(unreliable_manifest, "test", envir = new.env()),
  "unreliable uncertainty"
)
stopifnot(!isTRUE(unreliable_info$legacy))

# Explicit legacy loading checks only the supported top-level object shape.
BurdenMLE_DN_models_autism_test <- list(list(unclass(raw_fit)))
legacy_path <- file.path(fixture_dir, "legacy.Rdata")
save(BurdenMLE_DN_models_autism_test, autism_data, file = legacy_path)
legacy_environment <- new.env()
legacy_info <- expect_warning(
  load_model_artifact(
    artifact = "test", envir = legacy_environment, legacy_path = legacy_path
  ),
  "unmanifested legacy"
)
stopifnot(
  isTRUE(legacy_info$legacy),
  exists("BurdenMLE_DN_models_autism_test", legacy_environment,
         inherits = FALSE)
)

# Active consumers use the shared loader and contain no lexical model search.
consumer_files <- c(
  "summarize_main_autism.R", "summarize_cohort_sex.R", "summarize_ddid.R",
  "summarize_secondary_autism.R", "summarize_posterior_diagnostics.R",
  "summarize_autism_dd.R", "summarize_autism_dd_prevalence.R",
  "summarize_no_overlap.R", "forecasting_script_revision.R",
  "describe_us_father_age.R", "make_supplementary_tables.R",
  "make_manuscript_estimate_table.R", "figure_modules/make_sensitivity_figures.R"
)
for (consumer_file in file.path("analysis", "scripts", consumer_files)) {
  source_text <- paste(readLines(consumer_file, warn = FALSE), collapse = "\n")
  stopifnot(
    grepl("load_model_artifact\\(", source_text),
    !grepl("list\\.files\\(", source_text),
    !grepl("models_.*\\\\.Rdata", source_text)
  )
  parse(consumer_file)
}

producer_text <- paste(
  readLines("analysis/scripts/run_models.R", warn = FALSE), collapse = "\n"
)
stopifnot(
  grepl("model_manifest_.*analysis_mode.*current_date", producer_text),
  grepl("write_model_manifest\\(", producer_text),
  all(vapply(c("autism", "ddd", "age", "no_ces", "no_overlap", "test"),
             function(key) grepl(
               paste0("artifact_paths\\[\"", key, "\"\\]"), producer_text
             ), logical(1)))
)

runner_lines <- readLines("analysis/reproduce_study.sh", warn = FALSE)
assert_manifest_route <- function(script, variable) {
  line <- grep(paste0("/", script, "\\\""), runner_lines)
  stopifnot(length(line) == 1L)
  routed_lines <- runner_lines[line:min(length(runner_lines), line + 3L)]
  stopifnot(any(grepl(variable, routed_lines, fixed = TRUE)))
}
main_consumers <- c(
  "summarize_main_autism.R", "summarize_cohort_sex.R",
  "forecasting_script_revision.R", "summarize_autism_dd.R",
  "summarize_autism_dd_prevalence.R", "summarize_secondary_autism.R",
  "summarize_posterior_diagnostics.R", "summarize_ddid.R",
  "describe_us_father_age.R", "make_supplementary_tables.R",
  "make_manuscript_estimate_table.R"
)
for (script in main_consumers) {
  assert_manifest_route(script, "MAIN_MODEL_MANIFEST")
}
assert_manifest_route("summarize_no_overlap.R", "NO_OVERLAP_MODEL_MANIFEST")
assert_manifest_route("make_supplementary_figures.R", "NO_CES_MODEL_MANIFEST")

cat("Model artifact manifest contract tests passed.\n")
