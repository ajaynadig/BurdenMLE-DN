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

estimate_mutvar_trio <- function(model,
                                 genetic_data,
                                 prevalence,
                                 gamma_scaling_factor = 1) {
  require_case_rate_for_mutvar(genetic_data)
  annotation_mutvar <- vapply(seq_len(ncol(model$features)), function(index) {
    mixing_weights <- model$delta[index, ]
    mixing_weights[mixing_weights < 0] <- 0

    endpoint_samples <- sample(
      model$component_endpoints,
      size = 100000,
      prob = mixing_weights,
      replace = TRUE
    )
    log_rr_samples <- runif(
      100000,
      min = pmin(0, endpoint_samples),
      max = pmax(0, endpoint_samples)
    ) * gamma_scaling_factor
    rr_samples <- exp(log_rr_samples)
    mutation_rate_samples <- sample(
      genetic_data$case_rate[model$features[, index] == 1],
      size = 100000,
      replace = TRUE
    )

    sum(model$features[, index]) *
      prevalence * mean((rr_samples - 1)^2 * 2 * mutation_rate_samples) /
      (1 - prevalence)
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
