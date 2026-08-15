library(BurdenMLEDN)

quiet_fit <- function(args) {
  fit <- NULL
  invisible(capture.output(fit <- do.call(BurdenMLE_DN, args)))
  fit
}

expect_error_matching <- function(expression, pattern) {
  result <- try(force(expression), silent = TRUE)
  stopifnot(
    inherits(result, "try-error"),
    grepl(pattern, as.character(result), fixed = TRUE)
  )
}

raw_data <- data.frame(
  case_count = c(0, 1, 2, 0, 3, 1, 4, 0),
  case_rate = c(0.01, 0.02, 0.03, 0.015, 0.025, 0.035, 0.04, 0.012),
  N = rep(10, 8),
  row.names = paste0("g", seq_len(8))
)
expected_count <- 2 * raw_data$N * raw_data$case_rate
expected_data <- data.frame(
  case_count = raw_data$case_count,
  expected_count = expected_count,
  row.names = rownames(raw_data)
)
features <- matrix(
  c(rep(1, 4), rep(0, 4), rep(0, 4), rep(1, 4)),
  nrow = 8,
  ncol = 2,
  dimnames = list(rownames(raw_data), c("low", "high"))
)
endpoints <- c(0, log(2), log(5), log(10))

base_args <- list(
  features = features,
  component_endpoints = endpoints,
  prevalence = 0.02,
  bootstrap = FALSE,
  null_sim = FALSE,
  mutvar_est = FALSE,
  estimate_posteriors = TRUE,
  estimate_effective_penetrance = FALSE
)

weight_references <- list(
  mixsqp = structure(c(
    1.91855486395595e-15, 1.38856075487069e-15,
    0.991921994732332, 1.45728815880909e-15,
    0.00807800526766401, 0.634371315323135,
    1.83356293844075e-15, 0.365628684676862
  ), dim = c(2L, 4L)),
  EM = structure(c(
    8.37904064100953e-105, 9.88131291682493e-324,
    0.989016361170575, 2.48082828847263e-271,
    0.0109836388294251, 0.63437131774482,
    1.0771418585496e-242, 0.36562868225518
  ), dim = c(2L, 4L))
)

for (optimizer in c("mixsqp", "EM")) {
  raw_args <- base_args
  raw_args$input_data <- raw_data
  raw_args$optimizer <- optimizer
  raw_args$max_iter <- 10000
  raw_args$tol <- 1e-10
  raw_fit <- quiet_fit(raw_args)

  expected_args <- raw_args
  expected_args$input_data <- expected_data
  expected_fit <- quiet_fit(expected_args)

  tolerance <- if (optimizer == "mixsqp") 1e-10 else 1e-12
  stopifnot(
    isTRUE(all.equal(
      unname(raw_fit$delta),
      weight_references[[optimizer]],
      tolerance = tolerance,
      check.attributes = FALSE
    )),
    isTRUE(all.equal(
      unname(expected_fit$delta),
      unname(raw_fit$delta),
      tolerance = tolerance,
      check.attributes = FALSE
    )),
    isTRUE(all.equal(expected_fit$ll, raw_fit$ll, tolerance = 1e-12)),
    isTRUE(all.equal(
      expected_fit$conditional_likelihood,
      raw_fit$conditional_likelihood,
      tolerance = 0
    )),
    isTRUE(all.equal(
      expected_fit$posterior_gene_estimates,
      raw_fit$posterior_gene_estimates,
      tolerance = 1e-12
    )),
    identical(
      expected_fit$posterior_gene_estimates$Expected_Count,
      expected_count
    ),
    identical(expected_fit$input_summary$sample_size, NA_real_),
    identical(expected_fit$input_summary$case_rate_available, FALSE),
    identical(raw_fit$input_summary$sample_size, 10),
    identical(raw_fit$input_summary$case_rate_available, TRUE)
  )
}

process_data <- BurdenMLEDN:::process_data_trio
processed_expected <- process_data(expected_data)
stopifnot(identical(processed_expected$expected_count, expected_count))

# Dual-schema inputs retain the supplied expected counts exactly when they are
# inside the explicit absolute-plus-relative consistency tolerance.
within_tolerance <- 1e-12 + 1e-8 * abs(expected_count)
dual_data <- raw_data
dual_data$expected_count <- expected_count + 0.5 * within_tolerance
processed_dual <- process_data(dual_data)
stopifnot(identical(processed_dual$expected_count, dual_data$expected_count))

outside_tolerance <- raw_data
outside_tolerance$expected_count <- expected_count + 2 * within_tolerance
expect_error_matching(
  process_data(outside_tolerance),
  "Supplied expected_count values are inconsistent"
)

# Supported schema and column-value failures are diagnosed before fitting.
invalid_inputs <- list(
  raw_data[, "case_count", drop = FALSE],
  raw_data[, c("case_count", "case_rate"), drop = FALSE],
  raw_data[, c("case_count", "N"), drop = FALSE],
  transform(expected_data, expected_count = as.character(expected_count)),
  transform(expected_data, expected_count = replace(expected_count, 1, NA_real_)),
  transform(expected_data, expected_count = replace(expected_count, 1, Inf)),
  transform(expected_data, expected_count = replace(expected_count, 1, -1)),
  transform(expected_data, case_rate = replace(raw_data$case_rate, 1, NA_real_)),
  transform(expected_data, case_rate = replace(raw_data$case_rate, 1, -1)),
  transform(expected_data, N = rep(c(10, 11), each = 4)),
  transform(expected_data, N = rep(0, nrow(expected_data)))
)
for (invalid_input in invalid_inputs) {
  stopifnot(inherits(try(process_data(invalid_input), silent = TRUE), "try-error"))
}

# The public mutational-variance gate runs before component validation or model
# initialization and gives caller-owned derivation guidance.
mutvar_args <- base_args
mutvar_args$input_data <- expected_data
mutvar_args$component_endpoints <- "not a numeric endpoint vector"
mutvar_args$mutvar_est <- TRUE
expect_error_matching(
  quiet_fit(mutvar_args),
  "Mutational variance requires case_rate"
)
mutvar_args$prevalence <- NULL
expect_error_matching(
  quiet_fit(mutvar_args),
  "Mutational variance requires case_rate"
)

# A supplied mutation rate supports mutational variance without requiring N;
# expected counts remain the likelihood exposure and sample size remains unknown.
expected_with_rate <- expected_data
expected_with_rate$case_rate <- raw_data$case_rate
rate_args <- base_args
rate_args$input_data <- expected_with_rate
rate_args$optimizer <- "EM"
rate_args$mutvar_est <- TRUE
rate_args$estimate_posteriors <- FALSE
rate_args$max_iter <- 10000
rate_args$tol <- 1e-10
set.seed(20260815)
rate_fit <- quiet_fit(rate_args)
stopifnot(
  is.finite(rate_fit$mutvar_output$total_mutvar),
  identical(rate_fit$input_summary$sample_size, NA_real_),
  identical(rate_fit$input_summary$case_rate_available, TRUE)
)

# N may accompany exact expected counts without case_rate and is represented
# honestly without enabling mutation-rate-dependent estimands.
expected_with_n <- expected_data
expected_with_n$N <- raw_data$N
n_args <- base_args
n_args$input_data <- expected_with_n
n_args$optimizer <- "EM"
n_fit <- quiet_fit(n_args)
stopifnot(
  identical(n_fit$input_summary$sample_size, 10),
  identical(n_fit$input_summary$case_rate_available, FALSE)
)

# Other expected-count-only paths remain available because they depend on the
# likelihood exposure, not on per-haploid mutation rates.
capability_args <- base_args
capability_args$input_data <- expected_data
capability_args$optimizer <- "EM"
capability_args$null_sim <- TRUE
capability_args$n_null <- 2
capability_args$estimate_effective_penetrance <- TRUE
capability_args$max_iter <- 10000
capability_args$tol <- 1e-10
set.seed(7)
capability_fit <- quiet_fit(capability_args)
stopifnot(
  length(capability_fit$null_delta) == 2,
  is.finite(capability_fit$penetrance$effective_penetrance)
)

# A governed-data-free gnomAD-shaped calibration fixture uses gene_id as the
# explicit package identity and preserves its supplied synonymous expectations.
calibration_data <- data.frame(
  case_count = c(10, 3, 8, 1),
  expected_count = c(8.25, 4.5, 7.75, 2.25),
  gene_id = paste0("ENSG", seq_len(4))
)
rownames(calibration_data) <- calibration_data$gene_id
calibration_args <- list(
  input_data = calibration_data,
  features = NULL,
  component_endpoints = seq(-2, 2, length.out = 9),
  mutvar_est = FALSE,
  bootstrap = FALSE,
  null_sim = FALSE,
  estimate_posteriors = TRUE,
  estimate_effective_penetrance = FALSE,
  optimizer = "mixsqp"
)
calibration_fit <- quiet_fit(calibration_args)
stopifnot(
  identical(rownames(calibration_fit$posterior_gene_estimates), calibration_data$gene_id),
  identical(
    calibration_fit$posterior_gene_estimates$Expected_Count,
    calibration_data$expected_count
  )
)

# Prevalence is conditionally required only when the active path uses it.
generated_endpoint_args <- calibration_args
generated_endpoint_args$component_endpoints <- NULL
expect_error_matching(
  quiet_fit(generated_endpoint_args),
  "prevalence must be one finite number"
)

mutvar_prevalence_args <- calibration_args
mutvar_prevalence_args$input_data <- expected_with_rate
mutvar_prevalence_args$features <- features
mutvar_prevalence_args$component_endpoints <- endpoints
mutvar_prevalence_args$mutvar_est <- TRUE
expect_error_matching(
  quiet_fit(mutvar_prevalence_args),
  "prevalence must be one finite number"
)

penetrance_prevalence_args <- calibration_args
penetrance_prevalence_args$input_data <- expected_data
penetrance_prevalence_args$features <- features
penetrance_prevalence_args$component_endpoints <- endpoints
penetrance_prevalence_args$estimate_effective_penetrance <- TRUE
expect_error_matching(
  quiet_fit(penetrance_prevalence_args),
  "prevalence must be one finite number"
)

cat("Expected-count input contract tests passed.\n")
