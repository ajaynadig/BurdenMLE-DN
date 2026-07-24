# Model construction and posterior summaries for BurdenMLE-DN.

choose_component_endpoints_trio <- function(component_endpoints,
                                            no_cpts,
                                            prevalence) {
  if (!is.null(component_endpoints)) return(component_endpoints)
  seq(0, log(1 / prevalence), length.out = no_cpts)
}

initialize_model <- function(likelihood_function,
                             genetic_data,
                             component_endpoints,
                             features,
                             grid_size) {
  conditional_likelihood <- likelihood_function(
    genetic_data, component_endpoints, grid_size
  )
  no_cpts <- length(component_endpoints)

  if (is.null(features)) {
    features <- matrix(1, nrow = nrow(genetic_data), ncol = 1)
    rownames(features) <- rownames(genetic_data)
  }

  list(
    component_endpoints = component_endpoints,
    delta = matrix(1 / no_cpts, nrow = ncol(features), ncol = no_cpts),
    conditional_likelihood = conditional_likelihood,
    features = features,
    grid_size = grid_size
  )
}

normalize_log_likelihood_rows <- function(log_likelihood) {
  row_max <- apply(log_likelihood, 1, max)
  if (any(!is.finite(row_max))) {
    stop("At least one gene has zero likelihood at every integration point.")
  }
  scaled_likelihood <- exp(sweep(log_likelihood, 1, row_max, "-"))
  scaled_likelihood / rowSums(scaled_likelihood)
}

posterior_expectation <- function(model,
                                  genetic_data,
                                  function_to_integrate,
                                  grid_size) {
  weights <- model$features %*% model$delta
  posteriors <- weights * model$conditional_likelihood
  posteriors <- posteriors / rowSums(posteriors)
  no_tests <- nrow(model$conditional_likelihood)
  no_cpts <- length(model$component_endpoints)
  conditional_expectations <- matrix(NA_real_, no_tests, no_cpts)
  mu_grid <- seq(0.05, 1, by = 1 / grid_size)
  case_count_matrix <- replicate(length(mu_grid), genetic_data$case_count)

  for (component in seq_len(no_cpts)) {
    rate <- genetic_data$expected_count * t(replicate(
      no_tests, exp(mu_grid * model$component_endpoints[component])
    ))
    grid_posteriors <- normalize_log_likelihood_rows(
      dpois(case_count_matrix, rate, log = TRUE)
    )
    function_values <- function_to_integrate(
      mu_grid * model$component_endpoints[component]
    )
    conditional_expectations[, component] <- drop(
      grid_posteriors %*% function_values
    )
  }

  rowSums(posteriors * conditional_expectations)
}

precompute_effective_penetrance_moments <- function(
    model,
    genetic_data,
    gamma_scaling_factor = 1) {
  mu_grid <- seq(0.05, 1, by = 1 / model$grid_size)
  no_genes <- nrow(genetic_data)
  no_cpts <- length(model$component_endpoints)
  numerator <- matrix(NA_real_, no_genes, no_cpts)
  denominator <- matrix(NA_real_, no_genes, no_cpts)
  case_count_matrix <- matrix(
    genetic_data$case_count, nrow = no_genes, ncol = length(mu_grid)
  )

  for (component in seq_len(no_cpts)) {
    effects <- mu_grid * model$component_endpoints[component]
    grid_posteriors <- normalize_log_likelihood_rows(
      dpois(
        case_count_matrix,
        genetic_data$expected_count %o% exp(effects),
        log = TRUE
      )
    )
    numerator[, component] <- drop(
      grid_posteriors %*% ((exp(gamma_scaling_factor * effects) - 1) *
        exp(gamma_scaling_factor * effects))
    )
    denominator[, component] <- drop(
      grid_posteriors %*% (exp(gamma_scaling_factor * effects) - 1)
    )
  }

  list(numerator = numerator, denominator = denominator)
}

effective_penetrance_from_moments <- function(model, moments, prevalence) {
  posteriors <- (model$features %*% model$delta) *
    model$conditional_likelihood
  posteriors <- posteriors / rowSums(posteriors)
  numerator <- mean(rowSums(posteriors * moments$numerator))
  denominator <- mean(rowSums(posteriors * moments$denominator))
  if (!is.finite(denominator) || denominator == 0) return(NA_real_)
  prevalence * numerator / denominator
}

effective_penetrance_func <- function(model,
                                      genetic_data,
                                      prevalence,
                                      gamma_scaling_factor = 1,
                                      moments = NULL) {
  if (is.null(moments)) {
    moments <- precompute_effective_penetrance_moments(
      model, genetic_data, gamma_scaling_factor
    )
  }
  effective_penetrance_from_moments(model, moments, prevalence)
}

posterior_gene_samples <- function(model, num_samples = 1) {
  posteriors <- (model$features %*% model$delta) *
    model$conditional_likelihood
  posteriors <- posteriors / rowSums(posteriors)
  samples <- sapply(seq_len(nrow(posteriors)), function(gene) {
    endpoints <- sample(
      model$component_endpoints,
      size = num_samples,
      prob = posteriors[gene, ],
      replace = TRUE
    )
    runif(num_samples, min = pmin(0, endpoints), max = pmax(0, endpoints))
  })
  t(samples)
}
