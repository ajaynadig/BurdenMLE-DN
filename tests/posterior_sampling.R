library(BurdenMLEDN)

initialize_model <- BurdenMLEDN:::initialize_model
poisson_uniform_likelihood <- BurdenMLEDN:::poisson_uniform_likelihood
component_posteriors <- BurdenMLEDN:::component_posterior_probabilities
conditioned_draw <- BurdenMLEDN:::sample_conditioned_log_rr

expect_error <- function(code) {
  result <- try(force(code), silent = TRUE)
  stopifnot(inherits(result, "try-error"))
}

make_fit <- function(data, endpoints, weights, grid_size = 100L) {
  fit <- initialize_model(
    poisson_uniform_likelihood,
    data,
    endpoints,
    features = NULL,
    grid_size = grid_size
  )
  fit$delta[,] <- weights
  fit$posterior_gene_estimates <- data.frame(
    Case_Count = data$case_count,
    Expected_Count = data$expected_count,
    row.names = rownames(data)
  )
  class(fit) <- c("BurdenMLEDN_fit", "list")
  fit
}

component_mean <- function(count, expected, endpoint, function_value) {
  if (endpoint == 0) return(function_value(0))
  lower <- min(0, endpoint)
  upper <- max(0, endpoint)
  mode <- if (count == 0) {
    lower
  } else {
    min(upper, max(lower, log(count) - log(expected)))
  }
  log_likelihood <- function(x) count * x - expected * exp(x)
  scaled_density <- function(x) exp(log_likelihood(x) - log_likelihood(mode))
  denominator <- integrate(
    scaled_density, lower, upper, rel.tol = 1e-12, subdivisions = 1000L
  )$value
  integrate(
    function(x) function_value(x) * scaled_density(x),
    lower,
    upper,
    rel.tol = 1e-12,
    subdivisions = 1000L
  )$value / denominator
}

posterior_mean <- function(fit, data, gene, function_value) {
  probabilities <- component_posteriors(fit)[gene, ]
  sum(probabilities * vapply(
    fit$component_endpoints,
    function(endpoint) component_mean(
      data$case_count[gene], data$expected_count[gene], endpoint,
      function_value
    ),
    numeric(1)
  ))
}

# The prepared object preserves fitted component probabilities and requires the
# exact fitted gene identities and order.
informative_data <- data.frame(
  case_count = c(8, 0),
  expected_count = c(1, 0.2),
  row.names = c("informative", "weak")
)
informative_fit <- make_fit(
  informative_data, c(0, log(20)), c(0, 1), grid_size = 200L
)
sampler <- posterior_gene_sampler(informative_fit, informative_data)
prepared_probabilities <- sampler$cumulative_component_probabilities
prepared_probabilities[, -1] <-
  sampler$cumulative_component_probabilities[, -1, drop = FALSE] -
  sampler$cumulative_component_probabilities[, -ncol(
    sampler$cumulative_component_probabilities
  ), drop = FALSE]
stopifnot(isTRUE(all.equal(
  unname(prepared_probabilities),
  unname(component_posteriors(informative_fit)),
  tolerance = 1e-14
)))

reordered <- informative_data[2:1, , drop = FALSE]
expect_error(posterior_gene_sampler(informative_fit, reordered))
renamed <- informative_data
rownames(renamed)[1] <- "other"
expect_error(posterior_gene_sampler(informative_fit, renamed))
altered_count <- informative_data
altered_count$case_count[1] <- 0
expect_error(posterior_gene_sampler(informative_fit, altered_count))
altered_expected <- informative_data
altered_expected$expected_count[1] <- 2
expect_error(posterior_gene_sampler(informative_fit, altered_expected))
unnamed_fit <- informative_fit
unnamed_fit$grid_size <- 3L
legacy_grid <- seq(0.05, 1, by = 1 / unnamed_fit$grid_size)
legacy_counts <- replicate(unnamed_fit$grid_size, informative_data$case_count)
unnamed_fit$conditional_likelihood <- vapply(
  unnamed_fit$component_endpoints,
  function(endpoint) rowMeans(dpois(
    legacy_counts,
    informative_data$expected_count %o% exp(legacy_grid * endpoint)
  )),
  numeric(nrow(informative_data))
)
unnamed_fit$conditional_log_likelihood <- NULL
unnamed_fit$likelihood_log_scale <- NULL
legacy_sampler <- posterior_gene_sampler(unclass(unnamed_fit), informative_data)
stopifnot(identical(legacy_sampler$gene_ids, rownames(informative_data)))
unidentified_fit <- unnamed_fit
rownames(unidentified_fit$features) <- NULL
expect_error(posterior_gene_sampler(unidentified_fit, informative_data))
expect_error(posterior_gene_sampler(list(), informative_data))

# Shape and RNG contracts never simplify one-sample results.
set.seed(931)
one_sample <- posterior_gene_samples(sampler)
stopifnot(
  identical(dim(one_sample), c(2L, 1L)),
  identical(rownames(one_sample), rownames(informative_data)),
  identical(colnames(one_sample), "sample_1")
)
set.seed(241)
first_draws <- posterior_gene_samples(sampler, 7L)
set.seed(241)
second_draws <- posterior_gene_samples(sampler, 7L)
stopifnot(
  identical(first_draws, second_draws),
  identical(colnames(first_draws), paste0("sample_", 1:7))
)
for (bad_count in list(0, -1, 1.5, NA_real_, Inf, "1", c(1, 2))) {
  expect_error(posterior_gene_samples(sampler, bad_count))
}
malformed <- sampler
malformed$component_endpoints <- malformed$component_endpoints[-1]
expect_error(posterior_gene_samples(malformed))

# The private vector kernel also retains dimensions and support for a single
# gene, positive and negative intervals, and an exact zero point mass.
set.seed(18)
kernel_draws <- conditioned_draw(
  case_count = c(2, 0, 3),
  expected_count = c(1, 2, 1),
  endpoint = c(log(5), -log(4), 0)
)
stopifnot(
  length(conditioned_draw(2, 1, log(5))) == 1L,
  kernel_draws[1] >= 0,
  kernel_draws[1] <= log(5),
  kernel_draws[2] >= -log(4),
  kernel_draws[2] <= 0,
  identical(kernel_draws[3], 0)
)

# The informative audit fixture agrees with numerical integration. The old
# uniform-within-component mean is materially wrong for this same target.
target_rr_mean <- posterior_mean(informative_fit, informative_data, 1, exp)
target_log_mean <- posterior_mean(informative_fit, informative_data, 1, identity)
target_log_second <- posterior_mean(
  informative_fit, informative_data, 1, function(x) x^2
)
target_threshold <- posterior_mean(
  informative_fit, informative_data, 1, function(x) as.numeric(exp(x) >= 5)
)
old_uniform_rr_mean <- expm1(log(20)) / log(20)
stopifnot(
  abs(target_rr_mean - 7.99) < 0.01,
  abs(old_uniform_rr_mean - target_rr_mean) > 1.5
)

set.seed(314159)
informative_draws <- posterior_gene_samples(sampler, 150000L)
strong_draws <- informative_draws[1, ]
weak_draws <- informative_draws[2, ]
weak_mean <- posterior_mean(informative_fit, informative_data, 2, identity)
weak_second <- posterior_mean(
  informative_fit, informative_data, 2, function(x) x^2
)
weak_threshold <- posterior_mean(
  informative_fit, informative_data, 2, function(x) as.numeric(exp(x) >= 5)
)
stopifnot(
  abs(mean(exp(strong_draws)) - target_rr_mean) < 0.04,
  abs(mean(strong_draws) - target_log_mean) < 0.01,
  abs(var(strong_draws) -
        (target_log_second - target_log_mean^2)) < 0.015,
  abs(mean(exp(strong_draws) >= 5) - target_threshold) < 0.006,
  abs(mean(weak_draws) - weak_mean) < 0.01,
  abs(var(weak_draws) - (weak_second - weak_mean^2)) < 0.015,
  abs(mean(exp(weak_draws) >= 5) - weak_threshold) < 0.006
)

# A mixed negative/zero/positive grid exercises component selection and all
# within-component supports against the same integration oracle.
mixed_data <- data.frame(
  case_count = c(0, 2),
  expected_count = c(2, 1),
  row.names = c("negative_signal", "positive_signal")
)
mixed_fit <- make_fit(
  mixed_data,
  c(-log(4), 0, log(10)),
  c(0.3, 0.2, 0.5),
  grid_size = 200L
)
mixed_sampler <- posterior_gene_sampler(mixed_fit, mixed_data)
set.seed(2718)
mixed_draws <- posterior_gene_samples(mixed_sampler, 150000L)
for (gene in seq_len(nrow(mixed_data))) {
  exact_mean <- posterior_mean(mixed_fit, mixed_data, gene, identity)
  exact_second <- posterior_mean(
    mixed_fit, mixed_data, gene, function(x) x^2
  )
  exact_threshold <- posterior_mean(
    mixed_fit, mixed_data, gene, function(x) as.numeric(x >= log(2))
  )
  stopifnot(
    abs(mean(mixed_draws[gene, ]) - exact_mean) < 0.012,
    abs(var(mixed_draws[gene, ]) - (exact_second - exact_mean^2)) < 0.02,
    abs(mean(mixed_draws[gene, ] >= log(2)) - exact_threshold) < 0.006
  )
}

point_fit <- mixed_fit
point_fit$delta[,] <- c(0, 1, 0)
point_sampler <- posterior_gene_sampler(point_fit, mixed_data)
stopifnot(all(posterior_gene_samples(point_sampler, 20L) == 0))

# A synthetic mini-forecast preserves task-level reproducibility while using
# only the package sampler for its latent gene effects.
mini_forecast <- function(task_index) {
  set.seed(1771L + task_index * 104729L)
  log_rr <- posterior_gene_samples(mixed_sampler, 1L)[, 1L]
  counts <- rpois(nrow(mixed_data), 0.5 * exp(log_rr))
  c(log_rr = unname(log_rr), counts = counts)
}
stopifnot(identical(mini_forecast(3L), mini_forecast(3L)))

forecast_path <- "analysis/scripts/forecasting_script_revision.R"
legacy_path <- "analysis/scripts/secondary_analysis_functions.R"
if (file.exists(forecast_path) && file.exists(legacy_path)) {
  forecast_source <- paste(readLines(forecast_path, warn = FALSE), collapse = "\n")
  legacy_source <- paste(readLines(legacy_path, warn = FALSE), collapse = "\n")
  stopifnot(
    grepl("posterior_gene_sampler\\(", forecast_source),
    grepl("posterior_gene_samples\\(", forecast_source),
    !grepl("posterior_sampler <- function", forecast_source, fixed = TRUE),
    !grepl("sample_log_rr <- function", forecast_source, fixed = TRUE),
    !grepl("BurdenMLE_DN_power_forecasting", legacy_source, fixed = TRUE)
  )

  forecast_expressions <- parse(forecast_path)
  forecast_environment <- new.env(parent = globalenv())
  for (function_name in c("forecast_once", "run_task")) {
    definition <- which(vapply(forecast_expressions, function(expression) {
      is.call(expression) && identical(expression[[1]], quote(`<-`)) &&
        identical(expression[[2]], as.name(function_name))
    }, logical(1)))
    stopifnot(length(definition) == 1L)
    eval(forecast_expressions[[definition]], envir = forecast_environment)
  }

  forecast_data <- mixed_data
  forecast_data$N <- rep(10, nrow(forecast_data))
  forecast_data$case_rate <- forecast_data$expected_count / (2 * forecast_data$N)
  forecast_environment$tasks <- data.frame(
    n_new_case = 5,
    iteration = 1L,
    scenario = "synthetic",
    stringsAsFactors = FALSE
  )
  forecast_environment$scenario_list <- list(
    synthetic = list(
      data = forecast_data,
      sampler = mixed_sampler,
      scaling = 0.8,
      group = "test"
    )
  )
  forecast_environment$seed <- 1771L
  actual_task <- forecast_environment$run_task(1L)

  set.seed(1771L + 104729L)
  expected_task <- forecast_environment$forecast_once(
    forecast_data, mixed_sampler, 5, 0.8
  )
  expected_task$scenario <- "synthetic"
  expected_task$scenario_group <- "test"
  expected_task$n_new_case <- 5
  expected_task$total_sample_size <- 15
  expected_task$iteration <- 1L
  stopifnot(identical(actual_task, expected_task))
}

# Forecast-scale preparation remains cached and a batched draw is comfortably
# practical. The generous bound is intended to catch algorithmic regressions,
# not machine-level timing variation.
no_genes <- 17395L
large_ids <- paste0("gene", seq_len(no_genes))
large_data <- data.frame(
  case_count = rep(0, no_genes),
  expected_count = rep(0.2, no_genes),
  row.names = large_ids
)
large_fit <- make_fit(
  large_data,
  endpoints = c(0, log(2), log(10)),
  weights = c(0.7, 0.2, 0.1),
  grid_size = 10L
)
large_sampler <- posterior_gene_sampler(large_fit, large_data)
elapsed <- system.time(large_samples <- posterior_gene_samples(
  large_sampler, num_samples = 2L
))[["elapsed"]]
stopifnot(
  identical(dim(large_samples), c(no_genes, 2L)),
  elapsed < 5
)

cat("Likelihood-conditioned posterior sampler tests passed.\n")
