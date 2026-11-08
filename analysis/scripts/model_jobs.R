# Model-level jobs preserve the existing within-model bootstrap implementation.
# Preparation draws samples once, before jobs are scheduled; workers draw none.
model_job_code <- function(repo_dir) {
  files <- c(sort(list.files(file.path(repo_dir, "R"), pattern = "\\.R$", full.names = TRUE), method = "radix"),
             file.path(repo_dir, "analysis/scripts/model_jobs.R"))
  stats::setNames(unname(tools::md5sum(files)),
                  sub(paste0(normalizePath(repo_dir), "/"), "", normalizePath(files), fixed = TRUE))
}

model_job_bootstraps <- function(ids, n, seed) {
  with_local_seed(seed, matrix(vapply(seq_len(n), function(i) {
    ids[sample(seq_along(ids), replace = TRUE)]
  }, character(length(ids))), nrow = length(ids), ncol = n))
}

model_setup_code <- function(repo_dir) {
  files <- c(list.files(file.path(repo_dir, "R"), pattern = "\\.R$", full.names = TRUE),
             list.files(file.path(repo_dir, "analysis/scripts"), pattern = "\\.R$", full.names = TRUE))
  files <- sort(files, method = "radix")
  stats::setNames(unname(tools::md5sum(files)), substring(files, nchar(repo_dir) + 2L))
}

read_model_setup <- function(path, config, repo_dir) {
  expected <- readLines(paste0(path, ".md5"), warn = FALSE)
  if (!identical(unname(tools::md5sum(path)), expected)) stop("Setup snapshot checksum mismatch.")
  snapshot <- readRDS(path)
  if (!identical(snapshot$config, config) ||
      !identical(snapshot$code, model_setup_code(repo_dir))) stop("Setup snapshot settings/code mismatch.")
  snapshot$objects
}

atomic_job_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(".partial-", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp)
  if (!file.rename(tmp, path)) stop("Cannot publish ", path)
}

prepare_model_job <- function(arguments, label, path, repo_dir) {
  specification <- list(arguments = arguments, label = label, code = model_job_code(repo_dir))
  if (file.exists(path)) {
    if (!identical(readRDS(path), specification)) {
      stop("Existing job differs in data/settings/code: ", path, ". Use a new job directory.")
    }
  } else atomic_job_rds(specification, path)
  invisible(path)
}

model_job_result <- function(path) sub("\\.rds$", "_result.rds", path)

read_model_job_result <- function(path, repo_dir) {
  specification <- readRDS(path)
  result_path <- model_job_result(path)
  if (!file.exists(result_path)) stop("Incomplete model job: ", path)
  result <- readRDS(result_path)
  if (!identical(result$specification_md5, unname(tools::md5sum(path))) ||
      !identical(specification$code, model_job_code(repo_dir))) {
    stop("Stale checkpoint for ", path)
  }
  result$fit
}

execute_model_job <- function(path, repo_dir) {
  path <- normalizePath(path, mustWork = TRUE)
  if (file.exists(model_job_result(path))) {
    read_model_job_result(path, repo_dir)
    message("Validated existing checkpoint: ", basename(path))
    return(invisible(model_job_result(path)))
  }
  lock <- paste0(path, ".lock")
  if (!dir.create(lock, showWarnings = FALSE)) stop("Job already locked: ", path)
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  specification <- readRDS(path)
  if (!identical(specification$code, model_job_code(repo_dir))) stop("Worker code differs from prepared job.")
  message("Fitting ", specification$label)
  fit <- do.call(BurdenMLE_DN, specification$arguments)
  if (!isTRUE(fit$fit_status$converged) || !isTRUE(fit$uncertainty_reliable)) {
    stop("Unreliable fit/uncertainty: ", specification$label)
  }
  atomic_job_rds(list(specification_md5 = unname(tools::md5sum(path)), fit = fit,
                      session = sessionInfo()), model_job_result(path))
  invisible(model_job_result(path))
}
