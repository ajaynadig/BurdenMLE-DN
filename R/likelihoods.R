# Likelihood calculations for BurdenMLE-DN inference.

midpoint_grid <- function(grid_size) {
  (seq_len(grid_size) - 0.5) / grid_size
}

row_log_sum_exp <- function(log_values) {
  log_values <- as.matrix(log_values)
  row_max <- apply(log_values, 1, max)
  result <- row_max
  finite_rows <- is.finite(row_max)

  if (any(finite_rows)) {
    centered <- sweep(
      log_values[finite_rows, , drop = FALSE],
      1,
      row_max[finite_rows],
      "-"
    )
    result[finite_rows] <- row_max[finite_rows] +
      log(rowSums(exp(centered)))
  }

  result
}

row_log_mean_exp <- function(log_values) {
  row_log_sum_exp(log_values) - log(ncol(log_values))
}

poisson_uniform_likelihood <- function(genetic_data,
                                       component_endpoints,
                                       grid_size) {
  no_cpts <- length(component_endpoints)
  no_genes <- nrow(genetic_data)
  mu_grid <- midpoint_grid(grid_size)
  case_count_matrix <- matrix(
    genetic_data$case_count,
    nrow = no_genes,
    ncol = grid_size
  )
  conditional_log_likelihood <- matrix(
    NA_real_,
    nrow = no_genes,
    ncol = no_cpts
  )

  for (component in seq_len(no_cpts)) {
    effects <- mu_grid * component_endpoints[component]
    rates <- genetic_data$expected_count %o% exp(effects)
    conditional_log_likelihood[, component] <- row_log_mean_exp(
      dpois(case_count_matrix, rates, log = TRUE)
    )
  }

  likelihood_log_scale <- apply(conditional_log_likelihood, 1, max)
  if (any(!is.finite(likelihood_log_scale))) {
    stop("At least one gene has zero likelihood under every component.")
  }

  conditional_likelihood <- exp(sweep(
    conditional_log_likelihood,
    1,
    likelihood_log_scale,
    "-"
  ))
  rownames(conditional_log_likelihood) <- rownames(genetic_data)
  rownames(conditional_likelihood) <- rownames(genetic_data)
  names(likelihood_log_scale) <- rownames(genetic_data)

  list(
    conditional_likelihood = conditional_likelihood,
    conditional_log_likelihood = conditional_log_likelihood,
    likelihood_log_scale = likelihood_log_scale
  )
}
