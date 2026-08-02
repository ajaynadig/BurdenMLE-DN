# MixSQP optimization for the original trio-model representation.
#
# These functions are deliberately separate from EM.R so the established EM
# implementation and its default behavior remain unchanged.

validate_mixsqp_features <- function(features) {
  if (!is.matrix(features)) {
    features <- as.matrix(features)
  }
  if (anyNA(features) || any(!is.finite(features))) {
    stop("MixSQP requires a finite, non-missing features matrix.")
  }
  if (!all(features %in% c(0, 1)) || !all(rowSums(features) == 1)) {
    stop(
      "MixSQP requires one-hot features: every row must contain exactly one 1 ",
      "and otherwise only 0s. Use optimizer = \"EM\" for other feature designs."
    )
  }
  invisible(features)
}

MixSQP_fit <- function(model,
                       control = list(),
                       return_likelihood = TRUE) {
  if (!requireNamespace("mixsqp", quietly = TRUE)) {
    stop(
      "optimizer = \"mixsqp\" requires the mixsqp package. Install it with ",
      "install.packages(\"mixsqp\")."
    )
  }

  features <- validate_mixsqp_features(model$features)
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
    fits[[stratum]] <- mixsqp::mixsqp(
      likelihood,
      x0 = x0,
      control = control
    )
    if (fits[[stratum]]$status != "converged to optimal solution") {
      stop(
        "MixSQP did not converge in feature stratum ", stratum,
        ": ", fits[[stratum]]$status
      )
    }
    delta[stratum, ] <- fits[[stratum]]$x
  }

  dimnames(delta) <- dimnames(model$delta)
  model$delta <- delta
  if (return_likelihood) {
    model$ll <- absolute_mixture_log_likelihood(model)
  }
  model$mixsqp_output <- fits
  model
}

bootstrap_MixSQP <- function(model,
                             n_boot,
                             bootstrap_samples = NULL,
                             bootstrap_seeds = NULL,
                             control = list()) {
  if (is.null(bootstrap_samples)) {
    if (is.null(bootstrap_seeds)) {
      bootstrap_seeds <- seq_len(n_boot)
    }
    bootstrap_samples <- sapply(bootstrap_seeds, function(seed) {
      set.seed(seed)
      sample(seq_len(nrow(model$conditional_likelihood)), replace = TRUE)
    })
    cat("...bootstrap MixSQP")
  } else {
    cat("...bootstrap MixSQP with user-specified samples")
  }

  bootstrap_delta <- pblapply(seq_len(n_boot), function(iter) {
    model_boot <- model
    sample_indices <- bootstrap_samples[, iter]
    model_boot <- subset_model_likelihood_rows(model_boot, sample_indices)
    model_boot$features <- model_boot$features[sample_indices, , drop = FALSE]
    # Begin every resampled fit from the same neutral, strictly positive state
    # as the full-data fit so that zero full-data weights can re-enter.
    model_boot$delta[,] <- 1 / ncol(model_boot$delta)
    MixSQP_fit(
      model_boot,
      control = control,
      return_likelihood = FALSE
    )$delta
  })

  list(
    bootstrap_delta = bootstrap_delta,
    bootstrap_samples = bootstrap_samples
  )
}

null_MixSQP_trio <- function(genetic_data,
                             model,
                             n_null,
                             grid_size,
                             control = list()) {
  cat("...null MixSQP")

  pblapply(seq_len(n_null), function(iter) {
    if (iter %% 20 == 0) {
      cat(paste0("...", iter))
    }
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
      return_likelihood = FALSE
    )$delta
  })
}
