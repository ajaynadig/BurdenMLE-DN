# Core Expectation-Maximization scripts
EM_fit <- function(model,
                   max_iter,
                   tol = 1e-6,
                   return_likelihood = TRUE,
                   allow_status_return = FALSE) {
  max_iter <- validate_positive_integer(max_iter, "max_iter")
  if (length(tol) != 1L || !is.numeric(tol) || is.na(tol) ||
      !is.finite(tol)) {
    stop("tol must be one finite number.")
  }

  delta_dimnames <- dimnames(model$delta)
  OLS_denom <- solve(t(model$features) %*% model$features)
  OLS_denom_t_features <- OLS_denom %*% t(model$features)
  likelihood_trace <- rep(NA_real_, max_iter)
  relative_change <- tol + 1

  weights <- model$features %*% model$delta
  posteriors <- weights * model$conditional_likelihood
  likelihood_trace[1L] <- absolute_mixture_log_likelihood(model)
  posteriors <- posteriors / rowSums(posteriors)
  model$delta <- OLS_denom_t_features %*% posteriors
  iteration <- 1L

  while (is.finite(likelihood_trace[iteration]) &&
         isTRUE(relative_change > tol)) {
    if (iteration >= max_iter) break

    weights <- model$features %*% model$delta
    posteriors <- weights * model$conditional_likelihood
    likelihood_trace[iteration + 1L] <- tryCatch(
      absolute_mixture_log_likelihood(model),
      error = function(error) NA_real_
    )
    relative_change <- abs(
      likelihood_trace[iteration + 1L] - likelihood_trace[iteration]
    ) / abs(likelihood_trace[iteration])

    posteriors <- posteriors / rowSums(posteriors)
    model$delta <- OLS_denom_t_features %*% posteriors
    iteration <- iteration + 1L
  }

  dimnames(model$delta) <- delta_dimnames
  final_likelihood <- tryCatch(
    absolute_mixture_log_likelihood(model),
    error = function(error) NA_real_
  )
  relative_change <- abs(
    final_likelihood - likelihood_trace[iteration]
  ) / abs(likelihood_trace[iteration])
  while (isTRUE(relative_change > tol) && iteration < max_iter) {
    likelihood_trace[iteration + 1L] <- final_likelihood
    weights <- model$features %*% model$delta
    posteriors <- weights * model$conditional_likelihood
    posteriors <- posteriors / rowSums(posteriors)
    model$delta <- OLS_denom_t_features %*% posteriors
    dimnames(model$delta) <- delta_dimnames
    iteration <- iteration + 1L
    final_likelihood <- tryCatch(
      absolute_mixture_log_likelihood(model),
      error = function(error) NA_real_
    )
    relative_change <- abs(
      final_likelihood - likelihood_trace[iteration]
    ) / abs(likelihood_trace[iteration])
  }
  converged <- is.finite(relative_change) && relative_change <= tol
  backend_status <- if (converged) {
    "converged"
  } else if (iteration >= max_iter) {
    "maximum iterations reached"
  } else {
    "nonfinite relative likelihood change"
  }
  model$fit_status <- new_fit_status(
    optimizer = "EM",
    weights = model$delta,
    log_likelihood = final_likelihood,
    converged = converged,
    backend_message = backend_status,
    iterations = iteration
  )
  model$em_output <- list(
    likelihood_trace = c(likelihood_trace[seq_len(iteration)], final_likelihood),
    relative_change = relative_change
  )
  if (return_likelihood) model$ll <- final_likelihood

  if (!allow_status_return && !model$fit_status$usable) {
    stop("EM returned an unusable fit (", model$fit_status$code, ").")
  }
  if (!allow_status_return && !model$fit_status$converged) {
    stop("EM did not converge: ", model$fit_status$backend_message, ".")
  }
  model
}

bootstrap_EM <- function(model,
                         bootstrap_indices,
                         replicate_ids,
                         max_iter,
                         tol = 0) {
  cat("...bootstrap EM")
  fit_status <- pblapply(seq_along(replicate_ids), function(iteration) {
    model_boot <- subset_model_likelihood_rows(
      model, bootstrap_indices[, iteration]
    )
    model_boot$features <- model$features[
      bootstrap_indices[, iteration], , drop = FALSE
    ]
    boot_output <- EM_fit(
      model = model_boot,
      max_iter = max_iter,
      tol = tol,
      return_likelihood = FALSE,
      allow_status_return = TRUE
    )
    if (iteration %% 20L == 0L) {
      cat(paste0("...", iteration, "...(",
                 boot_output$fit_status$iterations, " iters)..."))
    }
    boot_output$fit_status
  })
  names(fit_status) <- replicate_ids
  fit_status
}

null_EM_trio <- function(genetic_data,
                         model,
                         max_iter,
                         null_seeds,
                         replicate_ids,
                         grid_size,
                         tol) {
  cat("...null EM")
  fit_status <- pblapply(seq_along(replicate_ids), function(iteration) {
    if (iteration %% 20L == 0L) cat(paste0("...", iteration))
    with_local_seed(null_seeds[[iteration]], {
      genetic_data_null <- genetic_data
      genetic_data_null$case_count <- rpois(
        nrow(genetic_data), genetic_data$expected_count
      )
      model_null <- initialize_model(
        likelihood_function = poisson_uniform_likelihood,
        genetic_data = genetic_data_null,
        component_endpoints = model$component_endpoints,
        features = model$features,
        grid_size = grid_size
      )
      EM_fit(
        model_null,
        max_iter = max_iter,
        tol = tol,
        return_likelihood = FALSE,
        allow_status_return = TRUE
      )$fit_status
    })
  })
  names(fit_status) <- replicate_ids
  fit_status
}
