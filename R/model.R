# Model construction and posterior summaries for BurdenMLE-DN.

choose_component_endpoints_trio <- function(component_endpoints,
                                            no_cpts,
                                            prevalence) {
  if (!is.null(component_endpoints)) {
    if (!is.numeric(component_endpoints) || !is.null(dim(component_endpoints)) ||
        length(component_endpoints) < 2L || anyNA(component_endpoints) ||
        any(!is.finite(component_endpoints)) ||
        anyDuplicated(component_endpoints) || any(diff(component_endpoints) <= 0)) {
      stop(
        "component_endpoints must be a numeric vector of at least two finite, unique, strictly increasing values."
      )
    }
    return(component_endpoints)
  }
  if (length(no_cpts) != 1L || !is.numeric(no_cpts) || is.na(no_cpts) ||
      !is.finite(no_cpts) ||
      no_cpts < 2 || no_cpts != as.integer(no_cpts)) {
    stop("no_cpts must be one integer of at least 2.")
  }
  effect_size_grid(prevalence = prevalence, no_cpts = no_cpts)
}

initialize_model <- function(likelihood_function,
                             genetic_data,
                             component_endpoints,
                             features,
                             grid_size) {
  likelihood <- likelihood_function(
    genetic_data, component_endpoints, grid_size
  )
  if (is.list(likelihood)) {
    conditional_likelihood <- likelihood$conditional_likelihood
    conditional_log_likelihood <- likelihood$conditional_log_likelihood
    likelihood_log_scale <- likelihood$likelihood_log_scale
  } else {
    conditional_likelihood <- likelihood
    conditional_log_likelihood <- log(conditional_likelihood)
    likelihood_log_scale <- rep(0, nrow(conditional_likelihood))
  }
  no_cpts <- length(component_endpoints)

  if (is.null(features)) {
    features <- matrix(
      1,
      nrow = nrow(genetic_data),
      ncol = 1L,
      dimnames = list(rownames(genetic_data), "all_genes")
    )
  }

  list(
    component_endpoints = component_endpoints,
    delta = matrix(
      1 / no_cpts,
      nrow = ncol(features),
      ncol = no_cpts,
      dimnames = list(colnames(features), NULL)
    ),
    conditional_likelihood = conditional_likelihood,
    conditional_log_likelihood = conditional_log_likelihood,
    likelihood_log_scale = likelihood_log_scale,
    features = features,
    grid_size = grid_size
  )
}

normalize_log_likelihood_rows <- function(log_likelihood) {
  row_normalizer <- row_log_sum_exp(log_likelihood)
  if (any(!is.finite(row_normalizer))) {
    stop("At least one gene has zero likelihood at every integration point.")
  }
  exp(sweep(log_likelihood, 1, row_normalizer, "-"))
}

component_log_likelihood <- function(model) {
  if (!is.null(model$conditional_log_likelihood)) {
    log_likelihood <- model$conditional_log_likelihood
    target_names <- rownames(model$conditional_likelihood)
    log_names <- rownames(log_likelihood)
    if (!is.null(target_names) && !is.null(log_names) &&
        all(target_names %in% log_names)) {
      log_likelihood <- log_likelihood[
        match(target_names, log_names),
        ,
        drop = FALSE
      ]
    }
    if (identical(dim(log_likelihood), dim(model$conditional_likelihood))) {
      return(log_likelihood)
    }
  }

  log_likelihood <- log(model$conditional_likelihood)
  if (!is.null(model$likelihood_log_scale)) {
    likelihood_log_scale <- model$likelihood_log_scale
    target_names <- rownames(model$conditional_likelihood)
    scale_names <- names(likelihood_log_scale)
    if (!is.null(target_names) && !is.null(scale_names) &&
        all(target_names %in% scale_names)) {
      likelihood_log_scale <- likelihood_log_scale[
        match(target_names, scale_names)
      ]
    }
    if (length(likelihood_log_scale) != nrow(log_likelihood)) {
      likelihood_log_scale <- rep(0, nrow(log_likelihood))
    }
    log_likelihood <- sweep(
      log_likelihood,
      1,
      likelihood_log_scale,
      "+"
    )
  }
  log_likelihood
}

component_posterior_probabilities <- function(model) {
  weights <- model$features %*% model$delta
  if (anyNA(weights) || any(!is.finite(weights)) || any(weights < 0)) {
    stop("Mixture weights must be finite and nonnegative.")
  }

  log_weights <- matrix(-Inf, nrow = nrow(weights), ncol = ncol(weights))
  positive <- weights > 0
  log_weights[positive] <- log(weights[positive])
  log_joint <- log_weights + component_log_likelihood(model)
  normalize_log_likelihood_rows(log_joint)
}

absolute_mixture_log_likelihood <- function(model) {
  weights <- model$features %*% model$delta
  if (anyNA(weights) || any(!is.finite(weights)) || any(weights < 0)) {
    stop("Mixture weights must be finite and nonnegative.")
  }

  log_weights <- matrix(-Inf, nrow = nrow(weights), ncol = ncol(weights))
  positive <- weights > 0
  log_weights[positive] <- log(weights[positive])
  sum(row_log_sum_exp(log_weights + component_log_likelihood(model)))
}

subset_model_likelihood_rows <- function(model, rows) {
  model$conditional_likelihood <-
    model$conditional_likelihood[rows, , drop = FALSE]
  if (!is.null(model$conditional_log_likelihood)) {
    model$conditional_log_likelihood <-
      model$conditional_log_likelihood[rows, , drop = FALSE]
  }
  if (!is.null(model$likelihood_log_scale)) {
    model$likelihood_log_scale <- model$likelihood_log_scale[rows]
  }
  model
}

posterior_expectation <- function(model,
                                  genetic_data,
                                  function_to_integrate,
                                  grid_size) {
  posteriors <- component_posterior_probabilities(model)
  no_tests <- nrow(model$conditional_likelihood)
  no_cpts <- length(model$component_endpoints)
  conditional_expectations <- matrix(NA_real_, no_tests, no_cpts)
  mu_grid <- midpoint_grid(grid_size)
  case_count_matrix <- matrix(
    genetic_data$case_count,
    nrow = no_tests,
    ncol = length(mu_grid)
  )

  for (component in seq_len(no_cpts)) {
    rate <- genetic_data$expected_count %o%
      exp(mu_grid * model$component_endpoints[component])
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
  mu_grid <- midpoint_grid(model$grid_size)
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
  posteriors <- component_posterior_probabilities(model)
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
  posteriors <- component_posterior_probabilities(model)
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
