# The existing fit collections contain pooled, clinical and population scenarios.
# Select by scientific convention, never by a separate reporting-only multiplier.
analysis_prevalence <- function(data, model_name) {
  i <- match(model_name, data$loop_vars$names)
  if (is.na(i)) stop("Unknown model: ", model_name)
  convention <- if (grepl("GeneDx|GeneDX", model_name)) {
    "clinical"
  } else if (grepl("\\b(ASC|SPARK)\\b", model_name)) {
    "population"
  } else "pooled"
  factor_index <- match(convention, c("pooled", "clinical", "population"))
  if (length(data$prev_factors) != 3L) stop("Expected three ASD prevalence scenarios.")
  prevalence <- data$loop_vars$prevalences[i] * data$prev_factors[factor_index]
  if (!is.finite(prevalence) || prevalence <= 0 || prevalence >= 1) {
    stop("Invalid analysis prevalence for ", model_name)
  }
  list(convention = convention, factor_index = factor_index, prevalence = prevalence)
}

assert_fit_prevalence <- function(model, prevalence, label) {
  if (!isTRUE(all.equal(model$prevalence, prevalence, tolerance = 1e-12))) {
    stop("Saved-fit/reporting prevalence mismatch for ", label,
         ". Select the matching prevalence scenario; do not rescale cached penetrance.")
  }
  invisible(TRUE)
}
