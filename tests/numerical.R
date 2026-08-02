library(BurdenMLEDN)

midpoint_grid <- BurdenMLEDN:::midpoint_grid
poisson_uniform_likelihood <- BurdenMLEDN:::poisson_uniform_likelihood
initialize_model <- BurdenMLEDN:::initialize_model
component_posteriors <- BurdenMLEDN:::component_posterior_probabilities
posterior_expectation <- BurdenMLEDN:::posterior_expectation
absolute_ll <- BurdenMLEDN:::absolute_mixture_log_likelihood
subset_likelihood_rows <- BurdenMLEDN:::subset_model_likelihood_rows
MixSQP_fit <- BurdenMLEDN:::MixSQP_fit
EM_fit <- BurdenMLEDN:::EM_fit

assert_equal <- function(actual, expected, tolerance = 1e-12) {
  stopifnot(isTRUE(all.equal(
    actual,
    expected,
    tolerance = tolerance,
    check.attributes = FALSE
  )))
}

test_log_sum_exp <- function(x) {
  row_max <- apply(x, 1, max)
  row_max + log(rowSums(exp(sweep(x, 1, row_max, "-"))))
}

legacy_likelihood <- function(genetic_data, component_endpoints) {
  mu_grid <- seq(0.05, 1, by = 0.1)
  counts <- matrix(
    genetic_data$case_count,
    nrow = nrow(genetic_data),
    ncol = length(mu_grid)
  )
  vapply(component_endpoints, function(endpoint) {
    rates <- genetic_data$expected_count %o% exp(mu_grid * endpoint)
    rowMeans(dpois(counts, rates))
  }, numeric(nrow(genetic_data)))
}

# The canonical grid contains exactly the requested midpoint values.
for (grid_size in c(1L, 3L, 10L, 20L, 100L)) {
  grid <- midpoint_grid(grid_size)
  stopifnot(
    length(grid) == grid_size,
    all(grid > 0),
    all(grid < 1)
  )
  assert_equal(grid, (seq_len(grid_size) - 0.5) / grid_size, 0)
}
assert_equal(midpoint_grid(10), seq(0.05, 1, by = 0.1), 1e-15)

moderate_data <- data.frame(
  case_count = c(0, 1, 2, 3, 1, 4, 0, 2),
  expected_count = c(0.2, 0.4, 0.6, 0.8, 0.3, 1, 0.1, 0.5),
  row.names = paste0("gene", seq_len(8))
)
endpoints <- c(0, log(1.5), log(3), log(8))

# Exercise the complete likelihood and posterior path at every target grid
# size, promoting any recycling warning to an error and comparing against a
# test-local direct log-space quadrature calculation.
check_full_grid_path <- function(grid_size) {
  old_options <- options(warn = 2)
  on.exit(options(old_options), add = TRUE)

  grid <- (seq_len(grid_size) - 0.5) / grid_size
  counts <- matrix(
    moderate_data$case_count,
    nrow = nrow(moderate_data),
    ncol = grid_size
  )
  direct_component_log <- vapply(endpoints, function(endpoint) {
    rates <- moderate_data$expected_count %o% exp(grid * endpoint)
    test_log_sum_exp(dpois(counts, rates, log = TRUE)) - log(grid_size)
  }, numeric(nrow(moderate_data)))

  likelihood <- poisson_uniform_likelihood(
    moderate_data,
    endpoints,
    grid_size
  )
  assert_equal(
    likelihood$conditional_log_likelihood,
    direct_component_log,
    1e-12
  )

  full_model <- initialize_model(
    poisson_uniform_likelihood,
    moderate_data,
    endpoints,
    features = NULL,
    grid_size = grid_size
  )
  component_weights <- matrix(
    1 / length(endpoints),
    nrow = nrow(moderate_data),
    ncol = length(endpoints)
  )
  direct_log_joint <- log(component_weights) + direct_component_log
  direct_log_normalizer <- test_log_sum_exp(direct_log_joint)
  direct_component_posterior <- exp(sweep(
    direct_log_joint,
    1,
    direct_log_normalizer,
    "-"
  ))
  direct_conditional_mean <- vapply(endpoints, function(endpoint) {
    effects <- grid * endpoint
    grid_log_likelihood <- dpois(
      counts,
      moderate_data$expected_count %o% exp(effects),
      log = TRUE
    )
    grid_log_normalizer <- test_log_sum_exp(grid_log_likelihood)
    drop(exp(sweep(
      grid_log_likelihood,
      1,
      grid_log_normalizer,
      "-"
    )) %*% exp(effects))
  }, numeric(nrow(moderate_data)))

  assert_equal(
    posterior_expectation(full_model, moderate_data, exp, grid_size),
    rowSums(direct_component_posterior * direct_conditional_mean),
    1e-12
  )
  assert_equal(
    absolute_ll(full_model),
    sum(direct_log_normalizer),
    1e-12
  )
}

invisible(lapply(c(1L, 3L, 10L, 20L, 100L), check_full_grid_path))

# Log-space integration reproduces the historical ordinary-scale calculation
# on moderate default-grid inputs.
stable_likelihood <- poisson_uniform_likelihood(
  moderate_data,
  endpoints,
  grid_size = 10
)
legacy_conditional <- legacy_likelihood(moderate_data, endpoints)
assert_equal(
  exp(stable_likelihood$conditional_log_likelihood),
  legacy_conditional,
  2e-14
)
stopifnot(
  all(is.finite(stable_likelihood$conditional_likelihood)),
  all(stable_likelihood$conditional_likelihood >= 0),
  all(rowSums(stable_likelihood$conditional_likelihood) > 0)
)

model <- initialize_model(
  poisson_uniform_likelihood,
  moderate_data,
  endpoints,
  features = NULL,
  grid_size = 10
)
legacy_model <- model
legacy_model$conditional_likelihood <- legacy_conditional
rownames(legacy_model$conditional_likelihood) <- rownames(moderate_data)
legacy_model$conditional_log_likelihood <- NULL
legacy_model$likelihood_log_scale <- NULL

# Repeated and reordered bootstrap rows keep every likelihood representation
# aligned, including row names and absolute offsets.
sample_rows <- c(6L, 2L, 6L, 1L)
resampled_model <- subset_likelihood_rows(model, sample_rows)
resampled_model$features <- model$features[sample_rows, , drop = FALSE]
assert_equal(
  resampled_model$conditional_likelihood,
  model$conditional_likelihood[sample_rows, , drop = FALSE],
  0
)
assert_equal(
  resampled_model$conditional_log_likelihood,
  model$conditional_log_likelihood[sample_rows, , drop = FALSE],
  0
)
assert_equal(
  resampled_model$likelihood_log_scale,
  model$likelihood_log_scale[sample_rows],
  0
)
assert_equal(
  component_posteriors(resampled_model),
  component_posteriors(model)[sample_rows, , drop = FALSE],
  1e-14
)
assert_equal(
  absolute_ll(resampled_model),
  sum(test_log_sum_exp(
    resampled_model$conditional_log_likelihood +
      log(matrix(
        1 / length(endpoints),
        nrow = length(sample_rows),
        ncol = length(endpoints)
      ))
  )),
  1e-12
)

# Both optimizer adapters preserve fitted weights, posterior summaries, and
# absolute likelihoods under row scaling.
mixsqp_fit <- MixSQP_fit(model)
legacy_mixsqp_fit <- MixSQP_fit(legacy_model)
assert_equal(mixsqp_fit$delta, legacy_mixsqp_fit$delta, 1e-10)
assert_equal(
  component_posteriors(mixsqp_fit),
  component_posteriors(legacy_mixsqp_fit),
  1e-11
)
assert_equal(mixsqp_fit$ll, legacy_mixsqp_fit$ll, 1e-11)
assert_equal(
  posterior_expectation(mixsqp_fit, moderate_data, exp, 10),
  posterior_expectation(legacy_mixsqp_fit, moderate_data, exp, 10),
  1e-11
)

# A negative tolerance forces the historical EM loop to its fixed iteration
# limit, separating numerical invariance from convergence-policy work.
em_fit <- EM_fit(model, max_iter = 20, tol = -1)
legacy_em_fit <- EM_fit(legacy_model, max_iter = 20, tol = -1)
assert_equal(em_fit$delta, legacy_em_fit$delta, 1e-12)
assert_equal(
  component_posteriors(em_fit),
  component_posteriors(legacy_em_fit),
  1e-12
)
assert_equal(em_fit$ll, legacy_em_fit$ll, 1e-11)
assert_equal(
  posterior_expectation(em_fit, moderate_data, exp, 10),
  posterior_expectation(legacy_em_fit, moderate_data, exp, 10),
  1e-11
)

# Legacy prepared models without the new fields remain usable.
legacy_model$delta <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 1)
modern_fixed <- model
modern_fixed$delta <- legacy_model$delta
assert_equal(
  component_posteriors(legacy_model),
  component_posteriors(modern_fixed),
  1e-12
)
assert_equal(
  posterior_expectation(legacy_model, moderate_data, exp, 10),
  posterior_expectation(modern_fixed, moderate_data, exp, 10),
  1e-12
)
assert_equal(
  absolute_ll(legacy_model),
  sum(log(drop(legacy_conditional %*% drop(legacy_model$delta)))),
  1e-12
)

# Extreme but finite Poisson data no longer collapse to zero likelihood rows.
extreme_data <- data.frame(
  case_count = c(200, 220),
  expected_count = c(0.01, 0.02),
  row.names = c("extreme1", "extreme2")
)
extreme_endpoints <- c(0, log(2), log(10))
extreme_model <- initialize_model(
  poisson_uniform_likelihood,
  extreme_data,
  extreme_endpoints,
  features = NULL,
  grid_size = 100
)
extreme_model$delta <- matrix(c(0.2, 0.3, 0.5), nrow = 1)
extreme_posteriors <- component_posteriors(extreme_model)
stopifnot(
  all(is.finite(extreme_model$conditional_likelihood)),
  all(is.finite(extreme_model$conditional_log_likelihood)),
  all(is.finite(extreme_model$likelihood_log_scale)),
  all(is.finite(extreme_posteriors)),
  all(is.finite(posterior_expectation(
    extreme_model,
    extreme_data,
    exp,
    100
  )))
)
assert_equal(rowSums(extreme_posteriors), rep(1, 2), 1e-14)

extreme_log_weights <- matrix(
  log(drop(extreme_model$delta)),
  nrow = nrow(extreme_data),
  ncol = length(extreme_endpoints),
  byrow = TRUE
)
direct_extreme_ll <- sum(test_log_sum_exp(
  extreme_model$conditional_log_likelihood + extreme_log_weights
))
assert_equal(absolute_ll(extreme_model), direct_extreme_ll, 1e-12)
stopifnot(is.finite(MixSQP_fit(extreme_model)$ll))
extreme_em <- EM_fit(extreme_model, max_iter = 10, tol = -1)
stopifnot(all(is.finite(extreme_em$ll)))

# A component-invariant row log shift leaves fitted weights and posteriors
# unchanged while changing every absolute log likelihood by the known sum.
row_shift <- c(7.25, -3, 0.5, 2, -1.75, 4, -0.25, 1.5)
shifted_model <- model
shifted_model$conditional_log_likelihood <- sweep(
  shifted_model$conditional_log_likelihood,
  1,
  row_shift,
  "+"
)
shifted_model$likelihood_log_scale <-
  shifted_model$likelihood_log_scale + row_shift

shifted_mixsqp <- MixSQP_fit(shifted_model)
assert_equal(shifted_mixsqp$delta, mixsqp_fit$delta, 1e-12)
assert_equal(
  component_posteriors(shifted_mixsqp),
  component_posteriors(mixsqp_fit),
  1e-12
)
assert_equal(shifted_mixsqp$ll, mixsqp_fit$ll + sum(row_shift), 1e-11)

shifted_em <- EM_fit(shifted_model, max_iter = 20, tol = -1)
assert_equal(shifted_em$delta, em_fit$delta, 1e-12)
assert_equal(
  component_posteriors(shifted_em),
  component_posteriors(em_fit),
  1e-12
)
assert_equal(
  shifted_em$ll,
  em_fit$ll + sum(row_shift),
  1e-11
)

cat("Numerical likelihood and posterior regression tests passed.\n")
