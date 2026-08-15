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

validate_positive_integer <- function(value, name) {
  if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
      !is.finite(value) || value < 1 || value != as.integer(value)) {
    stop(name, " must be one positive integer.")
  }
  as.integer(value)
}

resampling_seeds <- function(seed, replicates, name, prefix) {
  seed <- validate_positive_integer(seed, name)
  if (as.double(seed) + as.double(replicates) - 1 > .Machine$integer.max) {
    stop(name, " is too large to derive one seed per replicate.")
  }
  seeds <- as.integer(as.double(seed) + seq_len(replicates) - 1)
  names(seeds) <- paste0(prefix, seq_len(replicates))
  seeds
}

with_local_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

prepare_bootstrap_samples <- function(gene_ids,
                                      n_boot,
                                      bootstrap_samples,
                                      bootstrap_seed) {
  n_boot <- validate_positive_integer(n_boot, "n_boot")
  no_genes <- length(gene_ids)
  replicate_ids <- paste0("bootstrap_", seq_len(n_boot))

  if (is.null(bootstrap_samples)) {
    replicate_seeds <- resampling_seeds(
      bootstrap_seed, n_boot, "bootstrap_seed", "bootstrap_"
    )
    bootstrap_indices <- vapply(replicate_seeds, function(seed) {
      with_local_seed(seed, sample.int(no_genes, replace = TRUE))
    }, integer(no_genes))
    bootstrap_samples <- matrix(
      gene_ids[bootstrap_indices], nrow = no_genes, ncol = n_boot
    )
    source <- "generated"
  } else {
    if (!is.matrix(bootstrap_samples) || !is.character(bootstrap_samples)) {
      stop(
        "bootstrap_samples must be a character matrix of sampled gene IDs; ",
        "numeric row positions are no longer accepted."
      )
    }
    if (!identical(dim(bootstrap_samples), c(no_genes, n_boot))) {
      stop(
        "bootstrap_samples must have exactly one row per input gene and ",
        "n_boot columns."
      )
    }
    if (anyNA(bootstrap_samples) || any(bootstrap_samples == "")) {
      stop("bootstrap_samples must contain non-missing, nonempty gene IDs.")
    }
    unknown <- unique(bootstrap_samples[!bootstrap_samples %in% gene_ids])
    if (length(unknown) > 0L) {
      stop(
        "bootstrap_samples contains gene IDs absent from input_data: ",
        paste(utils::head(unknown, 5L), collapse = ", "), "."
      )
    }
    supplied_ids <- colnames(bootstrap_samples)
    if (!is.null(supplied_ids)) {
      if (anyNA(supplied_ids) || any(supplied_ids == "") ||
          anyDuplicated(supplied_ids)) {
        stop(
          "bootstrap_samples column names must be unique, nonempty, and ",
          "non-missing when supplied."
        )
      }
      replicate_ids <- supplied_ids
    }
    bootstrap_indices <- matrix(
      match(bootstrap_samples, gene_ids), nrow = no_genes, ncol = n_boot
    )
    replicate_seeds <- rep(NA_integer_, n_boot)
    source <- "supplied"
  }

  dimnames(bootstrap_samples) <- list(NULL, replicate_ids)
  dimnames(bootstrap_indices) <- list(NULL, replicate_ids)
  names(replicate_seeds) <- replicate_ids
  list(
    samples = bootstrap_samples,
    indices = bootstrap_indices,
    replicate_ids = replicate_ids,
    seeds = replicate_seeds,
    source = source
  )
}

mixture_weights_are_valid <- function(weights, tolerance = 1e-8) {
  is.matrix(weights) && is.numeric(weights) &&
    !anyNA(weights) && all(is.finite(weights)) &&
    all(weights >= 0) &&
    all(abs(rowSums(weights) - 1) <= tolerance)
}

new_fit_status <- function(optimizer,
                           weights,
                           log_likelihood,
                           converged,
                           backend_message,
                           iterations = NA_integer_) {
  likelihood_valid <- length(log_likelihood) == 1L &&
    is.numeric(log_likelihood) && !is.na(log_likelihood) &&
    is.finite(log_likelihood)
  weights_valid <- mixture_weights_are_valid(weights)
  usable <- likelihood_valid && weights_valid
  code <- if (!weights_valid) {
    "invalid_weights"
  } else if (!likelihood_valid) {
    "nonfinite_likelihood"
  } else if (isTRUE(converged)) {
    "converged"
  } else {
    "nonconverged"
  }
  list(
    optimizer = optimizer,
    log_likelihood = as.numeric(log_likelihood)[1L],
    weights = weights,
    converged = isTRUE(converged),
    usable = usable,
    code = code,
    backend_message = backend_message,
    iterations = as.integer(iterations)[1L]
  )
}

assess_uncertainty_fits <- function(fit_status, stage) {
  unusable <- names(fit_status)[
    !vapply(fit_status, function(record) isTRUE(record$usable), logical(1))
  ]
  if (length(unusable) > 0L) {
    stop(
      stage, " produced unusable replicate fits: ",
      paste(unusable, collapse = ", "), "."
    )
  }
  nonconverged <- names(fit_status)[
    !vapply(fit_status, function(record) isTRUE(record$converged), logical(1))
  ]
  if (length(nonconverged) > 0L) {
    warning(
      stage, " retained usable nonconverged replicate fits: ",
      paste(nonconverged, collapse = ", "), "."
    )
    return(FALSE)
  }
  TRUE
}

reconstruct_bootstrap_replicate <- function(model, genetic_data, iteration) {
  indices <- model$bootstrap_output$bootstrap_indices[, iteration]
  model_boot <- subset_model_likelihood_rows(model, indices)
  model_boot$features <- model$features[indices, , drop = FALSE]
  model_boot$delta <- model$bootstrap_output$bootstrap_delta[[iteration]]
  genetic_data_boot <- genetic_data[indices, , drop = FALSE]
  list(
    model = model_boot,
    genetic_data = genetic_data_boot,
    indices = indices
  )
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

#' Prepare a likelihood-conditioned gene posterior sampler
#'
#' Precomputes fitted component posterior probabilities for repeated sampling
#' of continuous gene log rate ratios. The fitted model is not modified or
#' refitted.
#'
#' @param fit A fitted `BurdenMLEDN_fit` object.
#' @param input_data The same gene-level data used to fit `fit`, with rows in
#'   exactly the same order. It must satisfy the input contract documented for
#'   [BurdenMLE_DN()].
#'
#' @return A prepared sampler for use with [posterior_gene_samples()].
#' @export
posterior_gene_sampler <- function(fit, input_data) {
  genetic_data <- process_data_trio(input_data)
  required_fields <- c(
    "component_endpoints", "delta", "features", "conditional_likelihood"
  )
  if (!is.list(fit) || any(!required_fields %in% names(fit))) {
    stop("fit must be a fitted BurdenMLEDN model.")
  }

  fit_gene_ids <- rownames(fit$features)
  if (is.null(fit_gene_ids) || anyNA(fit_gene_ids) ||
      any(!nzchar(fit_gene_ids)) || anyDuplicated(fit_gene_ids)) {
    stop("fit must contain unique, nonempty gene identities.")
  }
  likelihood_gene_ids <- rownames(fit$conditional_likelihood)
  if (!is.null(likelihood_gene_ids) &&
      !identical(likelihood_gene_ids, fit_gene_ids)) {
    stop("fit contains inconsistent gene identities.")
  }
  if (!identical(rownames(genetic_data), fit_gene_ids)) {
    stop(
      "input_data gene row names must be identical to the fitted genes in the same order."
    )
  }
  if (any(genetic_data$case_count > 0 & genetic_data$expected_count == 0)) {
    stop("Genes with positive observed counts must have positive expected counts.")
  }

  posteriors <- component_posterior_probabilities(fit)
  endpoints <- fit$component_endpoints
  if (!is.numeric(endpoints) || length(endpoints) != ncol(posteriors) ||
      anyNA(endpoints) || any(!is.finite(endpoints))) {
    stop("fit contains invalid component endpoints.")
  }
  if (nrow(posteriors) != nrow(genetic_data)) {
    stop("fit and input_data contain different numbers of genes.")
  }

  stored_data <- fit$posterior_gene_estimates
  if (is.data.frame(stored_data) &&
      all(c("Case_Count", "Expected_Count") %in% names(stored_data))) {
    stored_matches <- identical(rownames(stored_data), fit_gene_ids) &&
      identical(
        as.numeric(stored_data$Case_Count),
        as.numeric(genetic_data$case_count)
      ) &&
      isTRUE(all.equal(
        as.numeric(stored_data$Expected_Count),
        as.numeric(genetic_data$expected_count),
        tolerance = 1e-12
      ))
    if (!stored_matches) {
      stop("input_data counts do not match the data represented by fit.")
    }
  }

  cumulative <- t(apply(posteriors, 1L, cumsum))
  cumulative[cumulative < 0] <- 0
  cumulative[cumulative > 1] <- 1
  cumulative[, ncol(cumulative)] <- 1
  dimnames(cumulative) <- list(fit_gene_ids, NULL)

  structure(
    list(
      cumulative_component_probabilities = cumulative,
      case_count = genetic_data$case_count,
      expected_count = genetic_data$expected_count,
      component_endpoints = endpoints,
      gene_ids = fit_gene_ids
    ),
    class = c("BurdenMLEDN_posterior_sampler", "list")
  )
}

sample_conditioned_log_rr <- function(case_count, expected_count, endpoint) {
  if (!identical(length(case_count), length(expected_count)) ||
      !identical(length(case_count), length(endpoint))) {
    stop("Conditioned sampler inputs must have equal lengths.")
  }
  output <- numeric(length(endpoint))
  active <- endpoint != 0
  if (!any(active)) return(output)

  count <- case_count[active]
  expected <- expected_count[active]
  upper_endpoint <- pmax(0, endpoint[active])
  lower_endpoint <- pmin(0, endpoint[active])
  if (any(count > 0 & expected == 0)) {
    stop("Positive observed counts require positive expected counts.")
  }

  mode <- lower_endpoint
  positive_count <- count > 0
  mode[positive_count] <- pmin(
    upper_endpoint[positive_count],
    pmax(
      lower_endpoint[positive_count],
      log(count[positive_count]) - log(expected[positive_count])
    )
  )
  log_density <- function(x) count * x - expected * exp(x)
  maximum <- log_density(mode)

  remaining <- seq_along(count)
  accepted <- numeric(length(count))
  while (length(remaining)) {
    candidate <- runif(
      length(remaining),
      lower_endpoint[remaining],
      upper_endpoint[remaining]
    )
    log_ratio <- count[remaining] * candidate -
      expected[remaining] * exp(candidate) - maximum[remaining]
    keep <- log(runif(length(remaining))) <= pmin(0, log_ratio)
    if (any(keep)) accepted[remaining[keep]] <- candidate[keep]
    remaining <- remaining[!keep]
  }

  output[active] <- accepted
  output
}

#' Sample continuous gene effects from a fitted posterior
#'
#' Selects fitted mixture components from their gene-specific posterior
#' probabilities, then samples each continuous log-rate-ratio effect from the
#' Poisson likelihood conditioned on the selected component. Zero-width
#' components remain point masses at zero. The function uses R's ambient RNG
#' state and does not change or restore it.
#'
#' @param sampler A prepared sampler returned by [posterior_gene_sampler()].
#' @param num_samples Number of independent samples per gene.
#'
#' @return A numeric gene-by-sample matrix. Rows are named with fitted gene
#'   identities and columns are named `sample_1`, `sample_2`, and so on.
#' @export
posterior_gene_samples <- function(sampler, num_samples = 1) {
  if (!inherits(sampler, "BurdenMLEDN_posterior_sampler")) {
    stop("sampler must be created by posterior_gene_sampler().")
  }
  num_samples <- validate_positive_integer(num_samples, "num_samples")

  cumulative <- sampler$cumulative_component_probabilities
  no_genes <- length(sampler$gene_ids)
  valid <- is.matrix(cumulative) && is.numeric(cumulative) &&
    nrow(cumulative) == no_genes && ncol(cumulative) >= 1L &&
    !anyNA(cumulative) && all(is.finite(cumulative)) &&
    all(cumulative >= 0) && all(cumulative <= 1) &&
    all(cumulative[, ncol(cumulative)] == 1) &&
    (ncol(cumulative) == 1L || all(
      cumulative[, -1, drop = FALSE] >=
        cumulative[, -ncol(cumulative), drop = FALSE]
    )) &&
    is.character(sampler$gene_ids) && !anyNA(sampler$gene_ids) &&
    all(nzchar(sampler$gene_ids)) && !anyDuplicated(sampler$gene_ids) &&
    is.numeric(sampler$case_count) && is.numeric(sampler$expected_count) &&
    is.numeric(sampler$component_endpoints) &&
    length(sampler$case_count) == no_genes &&
    length(sampler$expected_count) == no_genes &&
    length(sampler$component_endpoints) == ncol(cumulative) &&
    !anyNA(sampler$case_count) && all(is.finite(sampler$case_count)) &&
    !anyNA(sampler$expected_count) && all(is.finite(sampler$expected_count)) &&
    !anyNA(sampler$component_endpoints) &&
    all(is.finite(sampler$component_endpoints))
  if (!valid) stop("sampler is malformed; recreate it with posterior_gene_sampler().")

  uniforms <- matrix(runif(no_genes * num_samples), nrow = no_genes)
  component <- vapply(seq_len(num_samples), function(sample_index) {
    as.integer(rowSums(cumulative <= uniforms[, sample_index]) + 1L)
  }, integer(no_genes))
  component <- matrix(component, nrow = no_genes, ncol = num_samples)
  endpoints <- matrix(
    sampler$component_endpoints[component],
    nrow = no_genes,
    ncol = num_samples
  )
  samples <- sample_conditioned_log_rr(
    rep(sampler$case_count, num_samples),
    rep(sampler$expected_count, num_samples),
    as.vector(endpoints)
  )
  samples <- matrix(samples, nrow = no_genes, ncol = num_samples)
  dimnames(samples) <- list(
    sampler$gene_ids,
    paste0("sample_", seq_len(num_samples))
  )
  samples
}
