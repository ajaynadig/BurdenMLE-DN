library(BurdenMLEDN)

quiet_fit <- function(args) {
  fit <- NULL
  invisible(capture.output(fit <- do.call(BurdenMLE_DN, args)))
  fit
}

expect_fit_error <- function(args) {
  result <- try(quiet_fit(args), silent = TRUE)
  stopifnot(inherits(result, "try-error"))
}

dat <- data.frame(
  case_count = c(0, 1, 2, 0, 3, 1, 4, 0),
  case_rate = c(0.01, 0.02, 0.03, 0.015, 0.025, 0.035, 0.04, 0.012),
  N = rep(10, 8),
  row.names = paste0("g", seq_len(8))
)
features <- matrix(
  c(rep(1, 4), rep(0, 4), rep(0, 4), rep(1, 4)),
  nrow = 8,
  ncol = 2,
  dimnames = list(rownames(dat), c("low", "high"))
)
endpoints <- c(0, log(2), log(5), log(10))

base_args <- list(
  input_data = dat,
  features = features,
  component_endpoints = endpoints,
  prevalence = 0.02,
  bootstrap = FALSE,
  null_sim = FALSE,
  mutvar_est = FALSE,
  estimate_posteriors = FALSE,
  estimate_effective_penetrance = FALSE
)

invalid_features <- list(
  continuous = { x <- features; x[1, ] <- c(0.5, 0.5); x },
  overlapping = { x <- features; x[1, ] <- 1; x },
  nonbinary = { x <- features; x[1, ] <- c(2, -1); x },
  all_zero_row = { x <- features; x[1, ] <- 0; x },
  empty_column = cbind(features, empty = 0),
  nonmatrix = as.data.frame(features),
  nonnumeric = { x <- matrix(as.character(features), 8, 2,
    dimnames = dimnames(features)); x },
  nonfinite = { x <- features; x[1, 1] <- Inf; x },
  unnamed_rows = { x <- features; rownames(x) <- NULL; x },
  duplicate_rows = { x <- features; rownames(x)[2] <- rownames(x)[1]; x },
  missing_rows = { x <- features; rownames(x)[2] <- NA_character_; x },
  reordered_rows = features[rev(seq_len(nrow(features))), , drop = FALSE],
  extra_gene = { x <- rbind(features, c(1, 0)); rownames(x)[9] <- "g9"; x },
  dimension_mismatch = features[-1, , drop = FALSE],
  unnamed_strata = { x <- features; colnames(x) <- NULL; x },
  empty_stratum_name = { x <- features; colnames(x)[1] <- ""; x },
  duplicate_strata = { x <- features; colnames(x) <- c("same", "same"); x },
  missing_stratum_name = { x <- features; colnames(x)[1] <- NA_character_; x }
)

invalid_inputs <- list(
  automatic = data.frame(
    case_count = dat$case_count,
    case_rate = dat$case_rate,
    N = dat$N
  ),
  missing = { x <- dat; attr(x, "row.names") <- c(NA_character_, rownames(dat)[-1]); x },
  empty = { x <- dat; attr(x, "row.names") <- c("", rownames(dat)[-1]); x },
  duplicate = { x <- dat; attr(x, "row.names") <- c("g1", "g1", rownames(dat)[-c(1, 2)]); x },
  one_gene = dat[1, , drop = FALSE]
)

for (optimizer in c("mixsqp", "EM")) {
  for (bad_features in invalid_features) {
    args <- base_args
    args$optimizer <- optimizer
    args$features <- bad_features
    expect_fit_error(args)
  }
  for (bad_input in invalid_inputs) {
    args <- base_args
    args$optimizer <- optimizer
    args$input_data <- bad_input
    args$features <- NULL
    expect_fit_error(args)
  }

  for (bad_endpoints in list(
    c("0", "1"),
    0,
    c(0, Inf),
    c(0, 2, 1),
    c(0, 1, 1)
  )) {
    args <- base_args
    args$optimizer <- optimizer
    args$component_endpoints <- bad_endpoints
    expect_fit_error(args)
  }

  for (bad_no_cpts in list(1, 2.5, NA_real_, "3", c(2, 3))) {
    args <- base_args
    args$optimizer <- optimizer
    args$component_endpoints <- NULL
    args$no_cpts <- bad_no_cpts
    expect_fit_error(args)
  }
}

# Explicit endpoints are the active specification and preserve their exact
# values and order even when no_cpts is invalid or missing.
precedence_args <- base_args
precedence_args$optimizer <- "EM"
precedence_args$no_cpts <- NA_real_
precedence_fit <- quiet_fit(precedence_args)
stopifnot(identical(precedence_fit$component_endpoints, endpoints))

missing_no_cpts_fit <- NULL
invisible(capture.output(
  missing_no_cpts_fit <- BurdenMLE_DN(
    dat,
    features = features,
    component_endpoints = endpoints,
    no_cpts =,
    prevalence = 0.02,
    bootstrap = FALSE,
    null_sim = FALSE,
    mutvar_est = FALSE,
    estimate_posteriors = FALSE,
    estimate_effective_penetrance = FALSE,
    optimizer = "EM"
  )
))
stopifnot(identical(missing_no_cpts_fit$component_endpoints, endpoints))

negative_endpoints <- c(-2, -0.5, 0.25)
negative_args <- base_args
negative_args$optimizer <- "EM"
negative_args$component_endpoints <- negative_endpoints
negative_fit <- quiet_fit(negative_args)
stopifnot(identical(negative_fit$component_endpoints, negative_endpoints))

# Frozen references from unchanged main before contract/name changes.
weight_references <- list(
  mixsqp = structure(c(
    1.91855486395595e-15, 1.38856075487069e-15,
    0.991921994732332, 1.45728815880909e-15,
    0.00807800526766401, 0.634371315323135,
    1.83356293844075e-15, 0.365628684676862
  ), dim = c(2L, 4L)),
  EM = structure(c(
    1.47749453399884e-06, 7.04600577078912e-25,
    0.888989988317672, 1.92467486032459e-17,
    0.111008534187794, 0.625124244673556,
    7.55552288009259e-17, 0.374875755326444
  ), dim = c(2L, 4L))
)

bootstrap_samples <- matrix(rep(seq_len(nrow(dat)), 2), ncol = 2)
processed_dat <- BurdenMLEDN:::process_data_trio(dat)
for (optimizer in c("mixsqp", "EM")) {
  args <- base_args
  args$optimizer <- optimizer
  args$bootstrap <- TRUE
  args$bootstrap_samples <- bootstrap_samples
  args$n_boot <- 2
  args$null_sim <- TRUE
  args$n_null <- 2
  args$mutvar_est <- TRUE
  args$max_iter <- 200
  args$max_iter_boot <- 200
  args$tol <- 1e-10
  set.seed(20260815)
  fit <- quiet_fit(args)

  stopifnot(isTRUE(all.equal(
    unname(fit$delta),
    weight_references[[optimizer]],
    tolerance = if (optimizer == "mixsqp") 1e-10 else 1e-12,
    check.attributes = FALSE
  )))
  gene_weights <- features %*% fit$delta
  stopifnot(
    all(is.finite(gene_weights)),
    all(gene_weights >= 0),
    isTRUE(all.equal(
      unname(rowSums(gene_weights)),
      rep(1, nrow(dat)),
      tolerance = 1e-8
    )),
    identical(rownames(fit$delta), colnames(features)),
    identical(rownames(mixture_weights(fit)), colnames(features)),
    identical(rownames(summary(fit)$mixture_weights), colnames(features)),
    identical(
      names(summary(fit)$mutational_variance$annot_mutvar),
      colnames(features)
    ),
    all(vapply(
      fit$bootstrap_output$bootstrap_delta,
      function(x) identical(rownames(x), colnames(features)),
      logical(1)
    )),
    all(vapply(
      fit$null_delta,
      function(x) identical(rownames(x), colnames(features)),
      logical(1)
    ))
  )

  named_vectors <- c("annot_mutvar", "frac_mutvar", "frac_expected", "enrichment")
  stopifnot(all(vapply(
    fit$mutvar_output[named_vectors],
    function(x) identical(names(x), colnames(features)),
    logical(1)
  )))
  named_intervals <- c("annot_mutvar_CI", "fracmutvar_CI", "enrich_CI")
  stopifnot(all(vapply(
    fit$mutvar_output[named_intervals],
    function(x) identical(colnames(x), colnames(features)),
    logical(1)
  )))

  unnamed_model <- fit
  colnames(unnamed_model$features) <- NULL
  rownames(unnamed_model$delta) <- NULL
  set.seed(99)
  named_mutvar <- BurdenMLEDN:::estimate_mutvar_trio(fit, processed_dat, 0.02)
  set.seed(99)
  unnamed_mutvar <- BurdenMLEDN:::estimate_mutvar_trio(
    unnamed_model, processed_dat, 0.02
  )
  stopifnot(isTRUE(all.equal(
    lapply(named_mutvar, unname),
    lapply(unnamed_mutvar, unname),
    tolerance = 0,
    check.attributes = FALSE
  )))
}

# The unstratified default has a stable public stratum identity.
unstratified_args <- base_args
unstratified_args$features <- NULL
unstratified_args$optimizer <- "EM"
unstratified_fit <- quiet_fit(unstratified_args)
stopifnot(
  identical(colnames(unstratified_fit$features), "all_genes"),
  identical(rownames(unstratified_fit$delta), "all_genes")
)

cat("Gene, feature, and component contract tests passed.\n")
