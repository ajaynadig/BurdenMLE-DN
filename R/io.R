# Input validation and preprocessing for de novo trio data.

process_data_trio <- function(input_data, features = NULL) {
  if (!is.data.frame(input_data)) {
    stop("input_data must be a data frame with one row per gene.")
  }
  if (!"case_count" %in% names(input_data)) {
    stop("input_data is missing required column: case_count.")
  }
  if (nrow(input_data) < 2L) {
    stop("input_data must contain at least two genes.")
  }

  gene_ids <- rownames(input_data)
  automatic_row_names <- .row_names_info(input_data, type = 1L) < 0L
  if (automatic_row_names || is.null(gene_ids) || anyNA(gene_ids) ||
      any(!nzchar(gene_ids)) || anyDuplicated(gene_ids)) {
    stop(
      "input_data must have explicit, unique, nonempty, non-missing gene row names."
    )
  }

  if (!is.numeric(input_data$case_count) || !is.null(dim(input_data$case_count)) ||
      anyNA(input_data$case_count) ||
      any(!is.finite(input_data$case_count)) ||
      any(input_data$case_count < 0) ||
      any(abs(input_data$case_count - round(input_data$case_count)) > 1e-8)) {
    stop("case_count must contain finite, non-missing, nonnegative integer counts.")
  }

  has_expected_count <- "expected_count" %in% names(input_data)
  has_case_rate <- "case_rate" %in% names(input_data)
  has_sample_size <- "N" %in% names(input_data)
  if (!has_expected_count && !(has_case_rate && has_sample_size)) {
    stop(
      "input_data must supply either (case_count, expected_count) or ",
      "(case_count, case_rate, N)."
    )
  }

  if (has_expected_count &&
      (!is.numeric(input_data$expected_count) ||
       !is.null(dim(input_data$expected_count)) ||
       anyNA(input_data$expected_count) ||
       any(!is.finite(input_data$expected_count)) ||
       any(input_data$expected_count < 0))) {
    stop("expected_count must contain finite, non-missing, nonnegative values.")
  }
  if (has_case_rate &&
      (!is.numeric(input_data$case_rate) || !is.null(dim(input_data$case_rate)) ||
       anyNA(input_data$case_rate) ||
       any(!is.finite(input_data$case_rate)) || any(input_data$case_rate < 0))) {
    stop(
      "case_rate must contain finite, non-missing, nonnegative per-haploid mutation rates."
    )
  }
  if (has_sample_size &&
      (!is.numeric(input_data$N) || !is.null(dim(input_data$N)) ||
       anyNA(input_data$N) ||
       any(!is.finite(input_data$N)) ||
       length(unique(input_data$N)) != 1L || any(input_data$N <= 0))) {
    stop("N must be one finite positive sample size repeated for every gene.")
  }

  if (has_case_rate && has_sample_size) {
    derived_expected_count <- 2 * input_data$N * input_data$case_rate
    if (any(!is.finite(derived_expected_count))) {
      stop("2 * N * case_rate must produce finite expected counts.")
    }
    if (has_expected_count) {
      consistency_tolerance <- 1e-12 + 1e-8 * pmax(
        abs(input_data$expected_count),
        abs(derived_expected_count)
      )
      if (any(abs(input_data$expected_count - derived_expected_count) >
              consistency_tolerance)) {
        stop(
          "Supplied expected_count values are inconsistent with ",
          "2 * N * case_rate."
        )
      }
    } else {
      input_data$expected_count <- derived_expected_count
    }
  }

  input_data
}

validate_features_trio <- function(input_data, features = NULL) {
  gene_ids <- rownames(input_data)
  if (is.null(features)) {
    features <- matrix(
      1,
      nrow = nrow(input_data),
      ncol = 1L,
      dimnames = list(gene_ids, "all_genes")
    )
    return(features)
  }

  if (!is.matrix(features) || !is.numeric(features)) {
    stop("features must be a numeric matrix.")
  }
  if (nrow(features) != nrow(input_data) || ncol(features) < 1L) {
    stop("features must have one row per input gene and at least one column.")
  }

  feature_gene_ids <- rownames(features)
  if (is.null(feature_gene_ids) || anyNA(feature_gene_ids) ||
      any(!nzchar(feature_gene_ids)) || anyDuplicated(feature_gene_ids)) {
    stop(
      "features must have unique, nonempty, non-missing gene row names."
    )
  }
  if (!identical(feature_gene_ids, gene_ids)) {
    stop(
      "features row names must be identical to input_data gene row names in the same order."
    )
  }

  stratum_names <- colnames(features)
  if (is.null(stratum_names) || anyNA(stratum_names) ||
      any(!nzchar(stratum_names)) || anyDuplicated(stratum_names)) {
    stop("features must have unique, nonempty, non-missing stratum column names.")
  }
  if (anyNA(features) || any(!is.finite(features))) {
    stop("features must contain only finite, non-missing values.")
  }
  if (!all(features %in% c(0, 1)) || !all(rowSums(features) == 1)) {
    stop(
      "features must be one-hot: every row must contain exactly one 1 and otherwise only 0s."
    )
  }
  if (any(colSums(features) == 0)) {
    stop("features must not contain empty strata.")
  }

  features
}
