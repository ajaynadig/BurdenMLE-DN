library(BurdenMLEDN)

fit_quietly <- function(args, warnings = FALSE) {
  fit <- NULL
  run <- function() invisible(capture.output(fit <<- do.call(BurdenMLE_DN, args)))
  if (warnings) {
    withCallingHandlers(run(), warning = function(w) invokeRestart("muffleWarning"))
  } else run()
  fit
}

expect_error <- function(args, text) {
  result <- try(fit_quietly(args, warnings = TRUE), silent = TRUE)
  stopifnot(inherits(result, "try-error"), grepl(text, result, fixed = TRUE))
}

dat <- data.frame(
  case_count = c(0, 1, 2, 0, 3, 1, 4, 0),
  case_rate = c(.01, .02, .03, .015, .025, .035, .04, .012),
  N = 10,
  row.names = paste0("g", 1:8)
)
features <- matrix(
  c(rep(1, 4), rep(0, 4), rep(0, 4), rep(1, 4)), 8, 2,
  dimnames = list(rownames(dat), c("low", "high"))
)
base <- list(
  input_data = dat,
  features = features,
  component_endpoints = c(0, log(2), log(5), log(10)),
  prevalence = .02,
  bootstrap = FALSE,
  null_sim = FALSE,
  mutvar_est = FALSE,
  estimate_posteriors = FALSE,
  estimate_effective_penetrance = FALSE,
  max_iter = 10000,
  tol = 1e-10
)

# The two optimizers expose the same final-fit contract.
for (optimizer in c("mixsqp", "EM")) {
  args <- base
  args$optimizer <- optimizer
  fit <- fit_quietly(args)
  status <- fit$fit_status
  stopifnot(
    identical(status$optimizer, optimizer),
    identical(status$weights, fit$delta),
    isTRUE(status$converged),
    isTRUE(status$usable),
    identical(status$code, "converged"),
    length(status$log_likelihood) == 1L,
    isTRUE(all.equal(
      status$log_likelihood,
      BurdenMLEDN:::absolute_mixture_log_likelihood(fit),
      tolerance = 1e-12
    )),
    identical(fit$ll, status$log_likelihood),
    identical(summary(fit)$fit_status, status),
    is.na(fit$uncertainty_reliable)
  )
  if (optimizer == "EM") {
    trace <- fit$em_output$likelihood_trace
    stopifnot(
      identical(tail(trace, 1L), status$log_likelihood),
      isTRUE(all.equal(
        fit$em_output$relative_change,
        abs(trace[length(trace)] - trace[length(trace) - 1L]) /
          abs(trace[length(trace) - 1L]),
        tolerance = 0
      ))
    )
  }
}

# A full-data fit may not silently return after exhausting its iteration limit.
args <- base
args$optimizer <- "EM"
args$max_iter <- 1L
expect_error(args, "EM did not converge")
args <- base
args$optimizer <- "mixsqp"
args$mixsqp_control <- list(maxiter.sqp = 1L)
expect_error(args, "MixSQP did not converge")

# Supplied samples are gene-identified, allow replacement, and survive reorder.
samples <- cbind(
  first = rownames(dat)[c(1, 1, 2, 3, 5, 6, 7, 8)],
  second = rev(rownames(dat))
)
for (optimizer in c("mixsqp", "EM")) {
  args <- base
  args$optimizer <- optimizer
  args$bootstrap <- TRUE
  args$bootstrap_samples <- samples
  args$n_boot <- 2L
  args$max_iter_boot <- 10000
  args$bootstrap_seed <- NA_real_ # inactive when samples are supplied
  fit <- fit_quietly(args)
  expected_indices <- matrix(
    match(samples, rownames(dat)), 8, 2,
    dimnames = list(NULL, colnames(samples))
  )
  stopifnot(
    identical(fit$bootstrap_output$bootstrap_samples, samples),
    identical(fit$bootstrap_output$bootstrap_indices, expected_indices),
    identical(fit$bootstrap_output$replicate_ids, colnames(samples)),
    identical(fit$bootstrap_output$sample_source, "supplied"),
    all(is.na(fit$bootstrap_output$replicate_seeds)),
    identical(names(fit$bootstrap_output$fit_status), colnames(samples)),
    identical(names(fit$bootstrap_output$bootstrap_delta), colnames(samples)),
    isTRUE(fit$uncertainty_reliable)
  )

  reordered_args <- args
  reordered_args$input_data <- dat[8:1, , drop = FALSE]
  reordered_args$features <- features[8:1, , drop = FALSE]
  reordered_fit <- fit_quietly(reordered_args)
  stopifnot(
    identical(reordered_fit$bootstrap_output$bootstrap_samples, samples),
    !identical(reordered_fit$bootstrap_output$bootstrap_indices,
               fit$bootstrap_output$bootstrap_indices),
    isTRUE(all.equal(
      lapply(reordered_fit$bootstrap_output$bootstrap_delta, unname),
      lapply(fit$bootstrap_output$bootstrap_delta, unname),
      tolerance = if (optimizer == "EM") 1e-10 else 1e-8
    ))
  )
}

# The shared validation rejects positional and malformed sample matrices.
bad_samples <- list(
  matrix(1:8, 8, 1),
  data.frame(sample = rownames(dat)),
  matrix(rownames(dat)[-1], 7, 1),
  { x <- matrix(rownames(dat), 8, 1); x[1] <- NA; x },
  { x <- matrix(rownames(dat), 8, 1); x[1] <- ""; x },
  { x <- matrix(rownames(dat), 8, 1); x[1] <- "unknown"; x },
  matrix(rep(rownames(dat), 2), 8, 2,
         dimnames = list(NULL, c("same", "same"))),
  matrix(rownames(dat), 8, 1, dimnames = list(NULL, "")),
  matrix(rownames(dat), 8, 1, dimnames = list(NULL, NA_character_))
)
for (optimizer in c("mixsqp", "EM")) {
  for (bad in bad_samples) {
    args <- base
    args$optimizer <- optimizer
    args$bootstrap <- TRUE
    args$bootstrap_samples <- bad
    args$n_boot <- ncol(as.matrix(bad))
    expect_error(args, "bootstrap_samples")
  }
}

# Bootstrap and null seeds replay independently and restore ambient RNG state.
seed_args <- base
seed_args$optimizer <- "mixsqp"
seed_args$features <- NULL
seed_args$bootstrap <- TRUE
seed_args$n_boot <- 3L
seed_args$bootstrap_seed <- 71L
set.seed(9001)
caller_seed <- .Random.seed
fit_one <- fit_quietly(seed_args)
stopifnot(identical(.Random.seed, caller_seed))
invisible(runif(5))
fit_two <- fit_quietly(seed_args)
stopifnot(
  identical(fit_one$bootstrap_output$bootstrap_samples,
            fit_two$bootstrap_output$bootstrap_samples),
  identical(fit_one$bootstrap_output$bootstrap_delta,
            fit_two$bootstrap_output$bootstrap_delta),
  identical(fit_one$bootstrap_output$replicate_seeds,
            c(bootstrap_1 = 71L, bootstrap_2 = 72L, bootstrap_3 = 73L))
)
seed_args$bootstrap_seed <- 72L
fit_three <- fit_quietly(seed_args)
stopifnot(!identical(fit_one$bootstrap_output$bootstrap_samples,
                     fit_three$bootstrap_output$bootstrap_samples))

null_args <- base
null_args$optimizer <- "mixsqp"
null_args$features <- NULL
null_args$null_sim <- TRUE
null_args$n_null <- 2L
null_args$null_seed <- 91L
set.seed(42)
caller_seed <- .Random.seed
null_one <- fit_quietly(null_args)
stopifnot(identical(.Random.seed, caller_seed))
invisible(runif(7))
null_two <- fit_quietly(null_args)
stopifnot(
  identical(null_one$null_delta, null_two$null_delta),
  identical(null_one$null_output$replicate_seeds,
            c(null_1 = 91L, null_2 = 92L)),
  identical(names(null_one$null_output$fit_status), c("null_1", "null_2"))
)

# A finite nonconverged replicate contributes, warns once, and marks CIs.
args <- base
args$optimizer <- "EM"
args$features <- NULL
args$bootstrap <- TRUE
args$bootstrap_samples <- matrix(
  rep("g1", 8), 8, 1,
  dimnames = list(NULL, "slow")
)
args$n_boot <- 1L
args$max_iter_boot <- 1L
args$mutvar_est <- TRUE
args$estimate_effective_penetrance <- TRUE
warnings <- character()
fit <- withCallingHandlers(fit_quietly(args), warning = function(w) {
  warnings <<- c(warnings, conditionMessage(w))
  invokeRestart("muffleWarning")
})
stopifnot(
  length(warnings) == 1L,
  grepl("slow", warnings, fixed = TRUE),
  identical(fit$uncertainty_reliable, FALSE),
  isTRUE(fit$bootstrap_output$fit_status$slow$usable),
  !isTRUE(fit$bootstrap_output$fit_status$slow$converged),
  all(is.finite(fit$mutvar_output$mutvar_CI)),
  all(is.finite(fit$penetrance$effective_penetrance_CI)),
  identical(dim(fit$mutvar_output$annot_mutvar_CI), c(2L, 1L)),
  identical(colnames(fit$mutvar_output$annot_mutvar_CI), "all_genes")
)

# MixSQP uses the same usable-nonconvergence policy at the adapter boundary.
genetic_data <- BurdenMLEDN:::process_data_trio(dat)
mix_model <- BurdenMLEDN:::initialize_model(
  BurdenMLEDN:::poisson_uniform_likelihood,
  genetic_data,
  base$component_endpoints,
  features = NULL,
  grid_size = 10
)
mix_model <- BurdenMLEDN:::MixSQP_fit(mix_model)
mix_status <- suppressWarnings(BurdenMLEDN:::bootstrap_MixSQP(
  mix_model, matrix(1:8, 8, 1), "slow",
  control = list(maxiter.sqp = 1L)
))
stopifnot(
  isTRUE(mix_status$slow$usable),
  !isTRUE(mix_status$slow$converged),
  identical(suppressWarnings(BurdenMLEDN:::assess_uncertainty_fits(
    mix_status, "Bootstrap"
  )), FALSE)
)

# Nonfinite likelihoods and invalid simplex weights are unusable and fatal.
bad_model <- BurdenMLEDN:::initialize_model(
  BurdenMLEDN:::poisson_uniform_likelihood,
  genetic_data,
  base$component_endpoints,
  features = NULL,
  grid_size = 10
)
bad_model$conditional_likelihood[,] <- 0
bad_model$conditional_log_likelihood[,] <- -Inf
bad_status <- BurdenMLEDN:::bootstrap_EM(
  bad_model, matrix(1:8, 8, 1), "bad", max_iter = 3L
)
bad_result <- try(BurdenMLEDN:::assess_uncertainty_fits(
  bad_status, "Bootstrap"
), silent = TRUE)
stopifnot(
  inherits(bad_result, "try-error"),
  grepl("bad", bad_result, fixed = TRUE)
)

valid_weights <- matrix(c(.5, .5), 1)
bad_records <- list(
  BurdenMLEDN:::new_fit_status("EM", valid_weights, NA_real_, FALSE, "test"),
  BurdenMLEDN:::new_fit_status(
    "EM", matrix(c(1.1, -.1), 1), -1, FALSE, "test"
  )
)
for (record in bad_records) {
  result <- try(BurdenMLEDN:::assess_uncertainty_fits(
    list(replicate_1 = record), "Test"
  ), silent = TRUE)
  stopifnot(inherits(result, "try-error"))
}

# Locally scoped generation also preserves the absence of an RNG state.
if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
  rm(".Random.seed", envir = .GlobalEnv)
}
invisible(fit_quietly(seed_args))
stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))

cat("Fit, resampling, and uncertainty contract tests passed.\n")
