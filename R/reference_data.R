#' Construct a BurdenMLE-DN effect-size grid
#'
#' By default, the largest rate ratio is `1 / prevalence`. A larger fixed
#' upper bound can be useful when prevalence may be misspecified or when very
#' large effects are scientifically plausible.
#'
#' @param prevalence Population prevalence on the 0--1 scale.
#' @param no_cpts Number of mixture components.
#' @param max_effect_size Optional maximum rate ratio. If `NULL`, uses
#'   `1 / prevalence`.
#'
#' @return A numeric vector of upper component endpoints on the log-rate-ratio
#'   scale.
#' @export
#'
#' @examples
#' effect_size_grid(prevalence = 0.02)
#' effect_size_grid(prevalence = 0.02, max_effect_size = 100)
effect_size_grid <- function(prevalence,
                             no_cpts = 10,
                             max_effect_size = NULL) {
  if (length(prevalence) != 1L || !is.finite(prevalence) ||
      prevalence <= 0 || prevalence >= 1) {
    stop("prevalence must be one finite number strictly between 0 and 1.")
  }
  if (length(no_cpts) != 1L || no_cpts < 2 ||
      no_cpts != as.integer(no_cpts)) {
    stop("no_cpts must be one integer of at least 2.")
  }
  if (is.null(max_effect_size)) max_effect_size <- 1 / prevalence
  if (length(max_effect_size) != 1L || !is.finite(max_effect_size) ||
      max_effect_size <= 1) {
    stop("max_effect_size must be one finite number greater than 1.")
  }
  seq(0, log(max_effect_size), length.out = no_cpts)
}

burdenmle_reference_strata <- function() {
  c(
    "LOEUF1_mu1", "LOEUF1_mu2",
    "LOEUF2", "LOEUF3", "LOEUF4", "LOEUF5"
  )
}

burdenmle_reference_rate_columns <- function() {
  base <- c(
    "mu_snp_PTV", "mu_snp_Mis2", "mu_snp_Mis1", "mu_snp_Mis0",
    "mu_snp_Syn"
  )
  c(base, paste0("corrected_", base))
}

validate_burdenmle_gene_reference <- function(reference) {
  character_columns <- c("gene_id", "gene_symbol", "analysis_stratum")
  numeric_columns <- c(
    "chromosome", "loeuf", "loeuf_quintile", "mutation_rate_correction",
    burdenmle_reference_rate_columns()
  )
  required <- c(
    "gene_id", "gene_symbol", "chromosome", "loeuf", "loeuf_quintile",
    "analysis_stratum", "mutation_rate_correction",
    burdenmle_reference_rate_columns()
  )
  if (!is.data.frame(reference) || !identical(names(reference), required)) {
    stop("The bundled gene reference has an invalid column schema.")
  }
  if (!all(vapply(reference[character_columns], is.character, logical(1))) ||
      !all(vapply(reference[numeric_columns], is.numeric, logical(1)))) {
    stop("The bundled gene reference has invalid column types.")
  }
  if (anyNA(reference) ||
      any(!nzchar(reference$gene_id)) ||
      any(!nzchar(reference$gene_symbol)) ||
      anyDuplicated(reference$gene_id)) {
    stop("The bundled gene reference has missing or duplicate gene identities.")
  }
  if (!identical(reference$gene_id, sort(reference$gene_id))) {
    stop("The bundled gene reference gene IDs are not in canonical order.")
  }
  if (any(!is.finite(as.matrix(reference[numeric_columns])))) {
    stop("The bundled gene reference contains nonfinite numeric values.")
  }
  if (any(reference$chromosome != as.integer(reference$chromosome)) ||
      any(!reference$chromosome %in% seq_len(22L))) {
    stop("The bundled gene reference contains invalid chromosomes.")
  }
  if (any(reference$loeuf < 0) ||
      any(reference$loeuf_quintile != as.integer(reference$loeuf_quintile)) ||
      any(!reference$loeuf_quintile %in% seq_len(5L))) {
    stop("The bundled gene reference contains invalid LOEUF values.")
  }

  allowed_strata <- burdenmle_reference_strata()
  if (!setequal(unique(reference$analysis_stratum), allowed_strata)) {
    stop("The bundled gene reference contains invalid analysis strata.")
  }
  expected_quintile <- c(1L, 1L, 2L, 3L, 4L, 5L)[
    match(reference$analysis_stratum, allowed_strata)
  ]
  if (!identical(as.integer(reference$loeuf_quintile), expected_quintile)) {
    stop("The bundled gene reference has inconsistent LOEUF strata.")
  }

  rate_columns <- burdenmle_reference_rate_columns()
  if (any(as.matrix(reference[c("mutation_rate_correction", rate_columns)]) < 0)) {
    stop("The bundled gene reference contains negative mutation rates.")
  }
  base_rates <- rate_columns[!startsWith(rate_columns, "corrected_")]
  for (column in base_rates) {
    corrected <- paste0("corrected_", column)
    if (!isTRUE(all.equal(
      reference[[corrected]],
      reference[[column]] * reference$mutation_rate_correction,
      tolerance = 1e-12
    ))) {
      stop("The bundled gene reference has inconsistent corrected mutation rates.")
    }
  }

  cumulative_rate <- rowSums(reference[paste0("corrected_", base_rates)])
  rate_threshold <- unname(stats::quantile(cumulative_rate, 0.8))
  expected_stratum <- ifelse(
    reference$loeuf_quintile == 1L,
    ifelse(cumulative_rate > rate_threshold, "LOEUF1_mu1", "LOEUF1_mu2"),
    paste0("LOEUF", reference$loeuf_quintile)
  )
  if (!identical(reference$analysis_stratum, expected_stratum)) {
    stop("The bundled gene reference has inconsistent mutation-rate strata.")
  }
  invisible(reference)
}

validate_burdenmle_loeuf_strata <- function(strata, reference) {
  required <- c(
    "analysis_stratum", "description", "genes", "loeuf_min", "loeuf_max"
  )
  if (!is.data.frame(strata) || !identical(names(strata), required)) {
    stop("The bundled LOEUF summary has an invalid column schema.")
  }
  allowed_strata <- burdenmle_reference_strata()
  if (!is.character(strata$analysis_stratum) ||
      !is.character(strata$description) ||
      !is.numeric(strata$genes) ||
      !is.numeric(strata$loeuf_min) ||
      !is.numeric(strata$loeuf_max) ||
      anyNA(strata) || any(!nzchar(strata$description)) ||
      !identical(strata$analysis_stratum, allowed_strata) ||
      anyDuplicated(strata$analysis_stratum) ||
      any(!is.finite(as.matrix(strata[c("genes", "loeuf_min", "loeuf_max")]))) ||
      any(strata$genes < 0) || any(strata$genes != as.integer(strata$genes)) ||
      any(strata$loeuf_min > strata$loeuf_max)) {
    stop("The bundled LOEUF summary contains invalid values.")
  }

  expected_genes <- vapply(allowed_strata, function(stratum) {
    sum(reference$analysis_stratum == stratum)
  }, integer(1))
  expected_min <- vapply(allowed_strata, function(stratum) {
    min(reference$loeuf[reference$analysis_stratum == stratum])
  }, numeric(1))
  expected_max <- vapply(allowed_strata, function(stratum) {
    max(reference$loeuf[reference$analysis_stratum == stratum])
  }, numeric(1))
  if (!identical(as.integer(strata$genes), unname(expected_genes)) ||
      !isTRUE(all.equal(strata$loeuf_min, unname(expected_min), tolerance = 1e-12)) ||
      !isTRUE(all.equal(strata$loeuf_max, unname(expected_max), tolerance = 1e-12))) {
    stop("The bundled LOEUF summary is not aligned with the gene reference.")
  }
  invisible(strata)
}

#' Load the BurdenMLE-DN gene reference
#'
#' Loads a public, count-free reference table containing gene identifiers,
#' LOEUF strata, and gene-level coding mutation rates used to construct the
#' manuscript's baseline annotation features. Logical schema, mutation-rate,
#' ordering, and stratum invariants are checked when the data are loaded. No
#' cohort variant or count data are included.
#'
#' @return A data frame with one row per gene.
#' @export
burdenmle_gene_reference <- function() {
  path <- system.file(
    "extdata", "burdenmle_gene_reference.csv", package = "BurdenMLEDN"
  )
  if (!nzchar(path)) stop("The installed gene reference could not be found.")
  reference <- utils::read.csv(
    path, stringsAsFactors = FALSE, check.names = FALSE
  )
  validate_burdenmle_gene_reference(reference)
  reference
}

#' Load the BurdenMLE-DN LOEUF strata summary
#'
#' @return A data frame describing the six baseline annotation strata, their
#'   gene counts, and observed LOEUF ranges in the bundled reference. Its
#'   counts and ranges are checked against the gene reference when loaded.
#' @export
burdenmle_loeuf_strata <- function() {
  path <- system.file(
    "extdata", "loeuf_strata_summary.csv", package = "BurdenMLEDN"
  )
  if (!nzchar(path)) stop("The installed LOEUF summary could not be found.")
  strata <- utils::read.csv(
    path, stringsAsFactors = FALSE, check.names = FALSE
  )
  reference <- burdenmle_gene_reference()
  validate_burdenmle_loeuf_strata(strata, reference)
  strata
}

#' Construct the manuscript LOEUF feature matrix
#'
#' Matches genes to the public BurdenMLE-DN reference and returns the six
#' one-hot annotation strata used in the main analyses: LOEUF quintiles 2--5
#' and the most constrained quintile split by cumulative mutation rate.
#'
#' @param gene_ids Character vector of Ensembl gene identifiers.
#' @param reference Reference data frame returned by
#'   [burdenmle_gene_reference()].
#'
#' @return A one-hot numeric matrix whose rows follow `gene_ids`.
#' @export
#'
#' @examples
#' reference <- burdenmle_gene_reference()
#' features <- loeuf_features(reference$gene_id[1:20], reference)
#' rowSums(features)
loeuf_features <- function(gene_ids,
                           reference = burdenmle_gene_reference()) {
  if (!is.character(gene_ids) || !length(gene_ids) ||
      anyNA(gene_ids) || anyDuplicated(gene_ids)) {
    stop("gene_ids must be a nonempty vector of unique, non-missing IDs.")
  }
  required <- c("gene_id", "analysis_stratum")
  missing_columns <- setdiff(required, names(reference))
  if (length(missing_columns)) {
    stop(
      "reference is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }
  index <- match(gene_ids, reference$gene_id)
  if (anyNA(index)) {
    stop(
      "The reference does not contain ", sum(is.na(index)),
      " requested gene ID(s)."
    )
  }
  levels <- burdenmle_reference_strata()
  stratum <- factor(reference$analysis_stratum[index], levels = levels)
  if (anyNA(stratum)) stop("The reference contains an unknown LOEUF stratum.")
  features <- stats::model.matrix(~ stratum - 1)
  colnames(features) <- levels
  rownames(features) <- gene_ids
  features
}
