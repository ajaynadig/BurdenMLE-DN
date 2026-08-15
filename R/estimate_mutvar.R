# Mutational-variance estimation for de novo burden models.

require_case_rate_for_mutvar <- function(genetic_data) {
  if (!"case_rate" %in% names(genetic_data)) {
    stop(
      "Mutational variance requires case_rate. For expected-count-only fitting, ",
      "set mutvar_est = FALSE. If N is known, derive case_rate = ",
      "expected_count / (2 * N) before fitting and rerun under your responsibility."
    )
  }
  invisible(TRUE)
}

uniform_squared_excess_rr_moment <- function(endpoint,
                                             gamma_scaling_factor = 1) {
  z <- gamma_scaling_factor * endpoint
  result <- numeric(length(z))
  nonzero <- z != 0
  if (any(nonzero)) {
    z_nonzero <- z[nonzero]
    result[nonzero] <-
      expm1(2 * z_nonzero) / (2 * z_nonzero) -
      2 * expm1(z_nonzero) / z_nonzero + 1
  }
  result
}

estimate_mutvar_trio <- function(model,
                                 genetic_data,
                                 prevalence,
                                 gamma_scaling_factor = 1) {
  require_case_rate_for_mutvar(genetic_data)
  component_moment <- uniform_squared_excess_rr_moment(
    model$component_endpoints,
    gamma_scaling_factor
  )

  annotation_mutvar <- vapply(seq_len(ncol(model$features)), function(index) {
    mixing_weights <- model$delta[index, ]
    mixing_weights[mixing_weights < 0] <- 0
    mixing_weights <- mixing_weights / sum(mixing_weights)
    rows <- model$features[, index] == 1

    prevalence / (1 - prevalence) *
      sum(2 * genetic_data$case_rate[rows]) *
      sum(mixing_weights * component_moment)
  }, numeric(1))
  names(annotation_mutvar) <- colnames(model$features)

  total_mutvar <- sum(annotation_mutvar)
  fraction_mutvar <- if (is.finite(total_mutvar) && total_mutvar > 0) {
    annotation_mutvar / total_mutvar
  } else {
    rep(NA_real_, length(annotation_mutvar))
  }
  names(fraction_mutvar) <- colnames(model$features)
  fraction_expected <- vapply(seq_len(ncol(model$features)), function(index) {
    sum(genetic_data$case_rate[model$features[, index] == 1])
  }, numeric(1)) / sum(genetic_data$case_rate)
  names(fraction_expected) <- colnames(model$features)

  list(
    total_mutvar = total_mutvar,
    annot_mutvar = annotation_mutvar,
    frac_mutvar = fraction_mutvar,
    frac_expected = fraction_expected,
    enrichment = fraction_mutvar / fraction_expected
  )
}
