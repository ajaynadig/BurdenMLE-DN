#' @export
print.BurdenMLEDN_fit <- function(x, ...) {
  cat("BurdenMLE-DN fit\n")
  cat("  Optimizer:", x$optimizer, "\n")
  if (!is.null(x$input_summary)) {
    cat("  Genes:", format(x$input_summary$genes, big.mark = ","), "\n")
    cat(
      "  Observed variants:",
      format(x$input_summary$observed_variants, big.mark = ","),
      "\n"
    )
  }
  if (!is.null(x$ll)) cat("  Log likelihood:", format(x$ll), "\n")
  if (!is.null(x$mutvar_output$total_mutvar)) {
    cat(
      "  Mutational variance:",
      format(x$mutvar_output$total_mutvar),
      "\n"
    )
  }
  invisible(x)
}

#' Summarize a BurdenMLE-DN fit
#'
#' @param object A fitted `BurdenMLEDN_fit` object.
#' @param ... Reserved for future use.
#' @return A compact list of model diagnostics and principal estimates.
#' @export
summary.BurdenMLEDN_fit <- function(object, ...) {
  list(
    optimizer = object$optimizer,
    input = object$input_summary,
    log_likelihood = object$ll,
    optimizer_elapsed = object$optimizer_elapsed,
    mixture_weights = object$delta,
    mutational_variance = object$mutvar_output,
    effective_penetrance = object$penetrance
  )
}

#' Extract fitted mixture weights
#'
#' @param object A fitted `BurdenMLEDN_fit` object.
#' @return A matrix of mixture weights, with one row per annotation stratum.
#' @export
mixture_weights <- function(object) {
  if (!inherits(object, "BurdenMLEDN_fit")) {
    stop("object must be a BurdenMLEDN_fit.")
  }
  object$delta
}
