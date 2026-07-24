# Input validation and preprocessing for de novo trio data.

process_data_trio <- function(input_data, features) {
  if (!is.data.frame(input_data)) {
    stop("input_data must be a data frame with one row per gene.")
  }
  required <- c("case_count", "case_rate", "N")
  missing_columns <- setdiff(required, names(input_data))
  if (length(missing_columns)) {
    stop(
      "input_data is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }
  if (!nrow(input_data)) stop("input_data contains no genes.")
  if (anyNA(input_data[, required, drop = FALSE]) ||
      any(!is.finite(as.matrix(input_data[, required, drop = FALSE])))) {
    stop("Required input columns must contain finite, non-missing values.")
  }
  if (any(input_data$case_count < 0) ||
      any(abs(input_data$case_count - round(input_data$case_count)) > 1e-8)) {
    stop("case_count must contain nonnegative integer counts.")
  }
  if (any(input_data$case_rate < 0)) {
    stop("case_rate must contain nonnegative per-haploid mutation rates.")
  }
  if (length(unique(input_data$N)) != 1L || any(input_data$N <= 0)) {
    stop("N must be one positive sample size repeated for every gene.")
  }
  if (is.null(rownames(input_data)) ||
      anyDuplicated(rownames(input_data))) {
    stop("input_data must have unique gene row names.")
  }
  input_data$expected_count <- 2 * input_data$N * input_data$case_rate

  if (!is.null(features) && any(is.na(features))) {
    stop("NAs present in features, please check")
  }
  if (!is.null(features) && nrow(features) != nrow(input_data)) {
    stop("features must have one row per input gene.")
  }
  if (!is.null(features) &&
      !all(rownames(input_data) == rownames(features))) {
    stop("features rownames do not match input data rownames, please check")
  }

  input_data
}
