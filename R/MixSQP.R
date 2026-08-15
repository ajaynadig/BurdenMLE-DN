# MixSQP optimization for the original trio-model representation.
#
# These functions are deliberately separate from EM.R so the established EM
# implementation and its default behavior remain unchanged.

MixSQP_fit <- function(model,
                       control = list(),
                       return_likelihood = TRUE,
                       allow_status_return = FALSE) {
  if (!requireNamespace("mixsqp", quietly = TRUE)) {
    stop(
      "optimizer = \"mixsqp\" requires the mixsqp package. Install it with ",
      "install.packages(\"mixsqp\")."
    )
  }

  features <- model$features
  conditional_likelihood <- model$conditional_likelihood
  no_cpts <- ncol(conditional_likelihood)

  # Begin optimization directly from the supplied positive weights. Keeping
  # the solution threshold at zero prevents numerically small components from
  # being irreversibly deleted before MixSQP has established whether they
  # belong in the optimum. Scale the active-set allowance with the user's
  # component grid, while retaining a floor suitable for the standard models.
  # An explicitly supplied control value still overrides this default.
  control <- utils::modifyList(
    list(
      verbose = FALSE,
      numiter.em = 0,
      maxiter.activeset = as.integer(max(20L, 2L * no_cpts)),
      zero.threshold.solution = 0
    ),
    control
  )
  if (anyNA(conditional_likelihood) || any(!is.finite(conditional_likelihood)) ||
      any(conditional_likelihood < 0)) {
    stop("MixSQP requires finite, nonnegative conditional likelihoods.")
  }

  no_strata <- ncol(features)
  delta <- matrix(NA_real_, nrow = no_strata, ncol = no_cpts)
  fits <- vector("list", no_strata)
  for (stratum in seq_len(no_strata)) {
    in_stratum <- features[, stratum] == 1
    if (!any(in_stratum)) {
      stop("Feature stratum ", stratum, " contains no genes.")
    }

    likelihood <- conditional_likelihood[in_stratum, , drop = FALSE]
    if (any(rowSums(likelihood) <= 0)) {
      stop(
        "Feature stratum ", stratum,
        " contains a gene with zero likelihood under every component."
      )
    }

    x0 <- model$delta[stratum, ]
    run_mixsqp <- function() mixsqp::mixsqp(
      likelihood, x0 = x0, control = control
    )
    fits[[stratum]] <- if (!allow_status_return) {
      run_mixsqp()
    } else {
      suppressWarnings(run_mixsqp())
    }
    delta[stratum, ] <- fits[[stratum]]$x
  }

  dimnames(delta) <- dimnames(model$delta)
  model$delta <- delta
  final_likelihood <- tryCatch(
    absolute_mixture_log_likelihood(model),
    error = function(error) NA_real_
  )
  converged_by_stratum <- vapply(
    fits,
    function(fit) identical(fit$status, "converged to optimal solution"),
    logical(1)
  )
  backend_status <- vapply(fits, function(fit) fit$status, character(1))
  names(backend_status) <- rownames(model$delta)
  iteration_count <- sum(vapply(fits, function(fit) {
    if (is.null(fit$progress)) 0L else nrow(fit$progress)
  }, integer(1)))
  model$fit_status <- new_fit_status(
    optimizer = "mixsqp",
    weights = model$delta,
    log_likelihood = final_likelihood,
    converged = all(converged_by_stratum),
    backend_message = paste(
      paste(names(backend_status), backend_status, sep = ": "),
      collapse = "; "
    ),
    iterations = iteration_count
  )
  if (return_likelihood) model$ll <- final_likelihood
  model$mixsqp_output <- fits
  if (!allow_status_return && !model$fit_status$usable) {
    stop("MixSQP returned an unusable fit (", model$fit_status$code, ").")
  }
  if (!allow_status_return && !model$fit_status$converged) {
    stop("MixSQP did not converge: ", model$fit_status$backend_message, ".")
  }
  model
}

bootstrap_MixSQP <- function(model,
                             bootstrap_indices,
                             replicate_ids,
                             control = list()) {
  cat("...bootstrap MixSQP")
  fit_status <- pblapply(seq_along(replicate_ids), function(iter) {
    model_boot <- model
    sample_indices <- bootstrap_indices[, iter]
    model_boot <- subset_model_likelihood_rows(model_boot, sample_indices)
    model_boot$features <- model_boot$features[sample_indices, , drop = FALSE]
    # Begin every resampled fit from the same neutral, strictly positive state
    # as the full-data fit so that zero full-data weights can re-enter.
    model_boot$delta[,] <- 1 / ncol(model_boot$delta)
    MixSQP_fit(
      model_boot,
      control = control,
      return_likelihood = FALSE,
      allow_status_return = TRUE
    )$fit_status
  })
  names(fit_status) <- replicate_ids
  fit_status
}

null_MixSQP_trio <- function(genetic_data,
                             model,
                             null_seeds,
                             replicate_ids,
                             grid_size,
                             control = list()) {
  cat("...null MixSQP")

  fit_status <- pblapply(seq_along(replicate_ids), function(iter) {
    if (iter %% 20 == 0) {
      cat(paste0("...", iter))
    }
    with_local_seed(null_seeds[[iter]], {
      genetic_data_null <- genetic_data
      genetic_data_null$case_count <- rpois(
        nrow(genetic_data),
        genetic_data$expected_count
      )
      model_null <- initialize_model(
        likelihood_function = poisson_uniform_likelihood,
        genetic_data = genetic_data_null,
        component_endpoints = model$component_endpoints,
        features = model$features,
        grid_size = grid_size
      )
      MixSQP_fit(
        model_null,
        control = control,
        return_likelihood = FALSE,
        allow_status_return = TRUE
      )$fit_status
    })
  })
  names(fit_status) <- replicate_ids
  fit_status
}
