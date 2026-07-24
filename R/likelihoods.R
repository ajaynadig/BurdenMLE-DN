# Likelihood calculations for BurdenMLE-DN inference.

poisson_uniform_likelihood <- function(genetic_data,
                                       component_endpoints,
                                       grid_size) {
  no_cpts <- length(component_endpoints)
  no_genes <- nrow(genetic_data)
  likelihood <- matrix(NA_real_, nrow = no_genes, ncol = no_cpts)
  mu_grid <- seq(0.05, 1, by = 1 / grid_size)
  case_count_matrix <- replicate(grid_size, genetic_data$case_count)

  for (component in seq_len(no_cpts)) {
    rate <- genetic_data$expected_count * t(replicate(
      no_genes,
      exp(mu_grid * component_endpoints[component])
    ))
    likelihood[, component] <- rowMeans(dpois(case_count_matrix, rate))
  }

  likelihood
}
