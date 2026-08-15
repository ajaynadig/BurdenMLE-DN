library(BurdenMLEDN)

moment <- BurdenMLEDN:::uniform_squared_excess_rr_moment
estimate_mutvar <- BurdenMLEDN:::estimate_mutvar_trio

# The analytic component moment agrees with direct integration for negative,
# zero, and positive endpoints, including gamma scaling.
endpoints <- c(-2, -0.25, 0, log(2), 3)
gamma <- 1.3
integrated_moments <- vapply(endpoints, function(endpoint) {
  integrate(
    function(unit_position) {
      expm1(gamma * endpoint * unit_position)^2
    },
    lower = 0,
    upper = 1,
    rel.tol = 1e-12
  )$value
}, numeric(1))
stopifnot(
  identical(moment(0), 0),
  isTRUE(all.equal(
    moment(endpoints, gamma),
    integrated_moments,
    tolerance = 1e-12
  ))
)

# A two-stratum fixture has hand-computable component and mutation-rate
# expectations. The reference formula is the exact implementation previously
# duplicated in the analysis scripts.
genes <- paste0("g", seq_len(6))
features <- matrix(
  c(rep(1, 3), rep(0, 3), rep(0, 3), rep(1, 3)),
  nrow = 6,
  ncol = 2,
  dimnames = list(genes, c("low", "high"))
)
delta <- matrix(
  c(0.2, -1e-14, 0.3, 0.4, 0.5, 0.6),
  nrow = 2,
  byrow = FALSE,
  dimnames = list(colnames(features), NULL)
)
component_endpoints <- c(-log(2), 0, log(3))
model <- list(
  features = features,
  delta = delta,
  component_endpoints = component_endpoints
)
genetic_data <- data.frame(
  case_count = c(0, 1, 2, 0, 1, 3),
  case_rate = c(0.01, 0.02, 0.04, 0.03, 0.05, 0.07),
  expected_count = c(0.2, 0.4, 0.8, 0.6, 1, 1.4),
  row.names = genes
)
prevalence <- 0.02
gamma <- 0.7

uniform_exp_mean <- function(endpoint, multiplier) {
  z <- multiplier * endpoint
  output <- rep(1, length(z))
  nonzero <- z != 0
  output[nonzero] <- expm1(z[nonzero]) / z[nonzero]
  output
}
reference_component_moment <-
  uniform_exp_mean(component_endpoints, 2 * gamma) -
  2 * uniform_exp_mean(component_endpoints, gamma) + 1
reference_annotation <- vapply(seq_len(ncol(features)), function(index) {
  rows <- features[, index] == 1
  weights <- pmax(delta[index, ], 0)
  weights <- weights / sum(weights)
  prevalence / (1 - prevalence) *
    sum(2 * genetic_data$case_rate[rows]) *
    sum(weights * reference_component_moment)
}, numeric(1))
names(reference_annotation) <- colnames(features)

exact <- estimate_mutvar(model, genetic_data, prevalence, gamma)
stopifnot(
  isTRUE(all.equal(
    exact$annot_mutvar,
    reference_annotation,
    tolerance = 1e-14
  )),
  isTRUE(all.equal(
    exact$total_mutvar,
    sum(reference_annotation),
    tolerance = 1e-14
  )),
  identical(names(exact$annot_mutvar), colnames(features)),
  identical(names(exact$frac_mutvar), colnames(features)),
  identical(names(exact$frac_expected), colnames(features)),
  identical(names(exact$enrichment), colnames(features))
)

one_stratum_model <- model
one_stratum_model$features <- matrix(
  1,
  nrow = nrow(features),
  ncol = 1,
  dimnames = list(rownames(features), "all_genes")
)
one_stratum_model$delta <- matrix(
  c(0.2, 0.3, 0.5),
  nrow = 1,
  dimnames = list("all_genes", NULL)
)
one_stratum <- estimate_mutvar(
  one_stratum_model,
  genetic_data,
  prevalence,
  gamma
)
stopifnot(
  length(one_stratum$annot_mutvar) == 1L,
  identical(names(one_stratum$annot_mutvar), "all_genes"),
  identical(
    estimate_mutvar(one_stratum_model, genetic_data, prevalence, 0)$total_mutvar,
    0
  )
)

# A high-volume simulation of the former estimator lands within its measured
# Monte Carlo sampling error of the exact result.
set.seed(20260815)
n_mc <- 500000L
stratum <- 1L
rows <- features[, stratum] == 1
weights <- pmax(delta[stratum, ], 0)
sampled_endpoint <- sample(
  component_endpoints,
  n_mc,
  replace = TRUE,
  prob = weights
)
sampled_log_rr <- runif(
  n_mc,
  pmin(0, sampled_endpoint),
  pmax(0, sampled_endpoint)
) * gamma
sampled_rate <- sample(genetic_data$case_rate[rows], n_mc, replace = TRUE)
sampled_contribution <- expm1(sampled_log_rr)^2 * 2 * sampled_rate
scale <- sum(rows) * prevalence / (1 - prevalence)
mc_estimate <- scale * mean(sampled_contribution)
mc_standard_error <- scale * sd(sampled_contribution) / sqrt(n_mc)
stopifnot(
  abs(mc_estimate - exact$annot_mutvar[stratum]) <= 6 * mc_standard_error
)

# Exact estimation is deterministic and does not create or advance RNG state.
if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
  rm(".Random.seed", envir = .GlobalEnv)
}
invisible(estimate_mutvar(model, genetic_data, prevalence, gamma))
stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))

set.seed(99)
rng_before <- .Random.seed
first <- estimate_mutvar(model, genetic_data, prevalence, gamma)
stopifnot(identical(.Random.seed, rng_before))
set.seed(12345)
second <- estimate_mutvar(model, genetic_data, prevalence, gamma)
stopifnot(identical(first, second))

# Point-mass components retain the existing zero-total convention.
zero_model <- model
zero_model$component_endpoints[] <- 0
zero <- estimate_mutvar(zero_model, genetic_data, prevalence)
stopifnot(
  identical(zero$total_mutvar, 0),
  identical(unname(zero$annot_mutvar), c(0, 0)),
  all(is.na(zero$frac_mutvar)),
  all(is.na(zero$enrichment)),
  identical(names(zero$annot_mutvar), colnames(features))
)

# The public full-fit and supplied-bootstrap paths are seed-independent, and
# replacing the derived estimand does not change the fitted mixture weights.
fit_data <- data.frame(
  case_count = c(0, 1, 2, 0, 3, 1, 4, 0),
  case_rate = c(0.01, 0.02, 0.03, 0.015, 0.025, 0.035, 0.04, 0.012),
  N = rep(10, 8),
  row.names = paste0("fit_g", seq_len(8))
)
fit_features <- matrix(
  c(rep(1, 4), rep(0, 4), rep(0, 4), rep(1, 4)),
  nrow = 8,
  ncol = 2,
  dimnames = list(rownames(fit_data), c("low", "high"))
)
bootstrap_samples <- cbind(
  seq_len(8),
  rev(seq_len(8)),
  c(1, 1, 2, 3, 5, 6, 7, 8)
)
fit_args <- list(
  input_data = fit_data,
  features = fit_features,
  component_endpoints = c(0, log(2), log(5), log(10)),
  prevalence = 0.02,
  bootstrap = TRUE,
  bootstrap_samples = bootstrap_samples,
  n_boot = ncol(bootstrap_samples),
  null_sim = FALSE,
  mutvar_est = TRUE,
  estimate_posteriors = FALSE,
  estimate_effective_penetrance = FALSE,
  optimizer = "EM",
  max_iter = 200,
  max_iter_boot = 200,
  tol = 1e-10
)
quiet_fit <- function(seed) {
  output <- NULL
  set.seed(seed)
  invisible(capture.output(output <- do.call(BurdenMLE_DN, fit_args)))
  output
}
fit_one <- quiet_fit(1)
fit_two <- quiet_fit(9876)
frozen_main_delta <- structure(c(
  1.47749453399884e-06, 7.04600577078912e-25,
  0.888989988317672, 1.92467486032459e-17,
  0.111008534187794, 0.625124244673556,
  7.55552288009259e-17, 0.374875755326444
), dim = c(2L, 4L))
stopifnot(
  isTRUE(all.equal(
    unname(fit_one$delta),
    frozen_main_delta,
    tolerance = 1e-12,
    check.attributes = FALSE
  )),
  identical(fit_one$delta, fit_two$delta),
  identical(fit_one$bootstrap_output$bootstrap_delta,
            fit_two$bootstrap_output$bootstrap_delta),
  identical(fit_one$mutvar_output, fit_two$mutvar_output),
  length(fit_one$mutvar_output$mutvar_CI) == 2L,
  identical(
    colnames(fit_one$mutvar_output$annot_mutvar_CI),
    colnames(fit_features)
  )
)

# Null datasets remain stochastic, but every fitted null weight matrix is
# passed through the same deterministic mutational-variance estimator.
null_args <- fit_args
null_args$bootstrap <- FALSE
null_args$bootstrap_samples <- NULL
null_args$n_boot <- 1
null_args$null_sim <- TRUE
null_args$n_null <- 3
set.seed(20260816)
null_fit <- NULL
invisible(capture.output(null_fit <- do.call(BurdenMLE_DN, null_args)))
processed_fit_data <- BurdenMLEDN:::process_data_trio(fit_data)
null_reference <- vapply(seq_along(null_fit$null_delta), function(index) {
  null_model <- null_fit
  null_model$delta <- null_fit$null_delta[[index]]
  estimate_mutvar(null_model, processed_fit_data, 0.02)$total_mutvar
}, numeric(1))
stopifnot(
  length(null_fit$mutvar_output$null_mutvar_ests) == null_args$n_null,
  identical(null_fit$mutvar_output$null_mutvar_ests, null_reference)
)

cat("Exact mutational-variance tests passed.\n")
