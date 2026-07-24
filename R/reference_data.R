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

#' Load the BurdenMLE-DN gene reference
#'
#' Loads a public, count-free reference table containing gene identifiers,
#' LOEUF strata, and gene-level coding mutation rates used to construct the
#' manuscript's baseline annotation features. No cohort variant or count data
#' are included.
#'
#' @return A data frame with one row per gene.
#' @export
burdenmle_gene_reference <- function() {
  path <- system.file(
    "extdata", "burdenmle_gene_reference.csv", package = "BurdenMLEDN"
  )
  if (!nzchar(path)) stop("The installed gene reference could not be found.")
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

#' Load the BurdenMLE-DN LOEUF strata summary
#'
#' @return A data frame describing the six baseline annotation strata, their
#'   gene counts, and observed LOEUF ranges in the bundled reference.
#' @export
burdenmle_loeuf_strata <- function() {
  path <- system.file(
    "extdata", "loeuf_strata_summary.csv", package = "BurdenMLEDN"
  )
  if (!nzchar(path)) stop("The installed LOEUF summary could not be found.")
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
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
  levels <- c(
    "LOEUF1_mu1", "LOEUF1_mu2",
    "LOEUF2", "LOEUF3", "LOEUF4", "LOEUF5"
  )
  stratum <- factor(reference$analysis_stratum[index], levels = levels)
  if (anyNA(stratum)) stop("The reference contains an unknown LOEUF stratum.")
  features <- stats::model.matrix(~ stratum - 1)
  colnames(features) <- levels
  rownames(features) <- gene_ids
  features
}
