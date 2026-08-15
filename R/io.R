# Input validation and preprocessing for de novo trio data.

process_data_trio <- function(input_data, features = NULL) {
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
  if (nrow(input_data) < 2L) {
    stop("input_data must contain at least two genes.")
  }
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
  gene_ids <- rownames(input_data)
  automatic_row_names <- .row_names_info(input_data, type = 1L) < 0L
  if (automatic_row_names || is.null(gene_ids) || anyNA(gene_ids) ||
      any(!nzchar(gene_ids)) || anyDuplicated(gene_ids)) {
    stop(
      "input_data must have explicit, unique, nonempty, non-missing gene row names."
    )
  }
  input_data$expected_count <- 2 * input_data$N * input_data$case_rate

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
