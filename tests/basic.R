library(BurdenMLEDN)

path <- system.file("extdata", "synthetic_example.csv", package = "BurdenMLEDN")
dat <- read.csv(path, row.names = "gene")

set.seed(1)
fit <- BurdenMLE_DN(
  dat,
  prevalence = 0.02,
  no_cpts = 6,
  bootstrap = FALSE,
  mutvar_est = FALSE,
  estimate_effective_penetrance = FALSE
)

stopifnot(
  inherits(fit, "BurdenMLEDN_fit"),
  identical(fit$optimizer, "mixsqp"),
  isTRUE(all.equal(sum(fit$delta), 1, tolerance = 1e-6)),
  is.finite(fit$ll),
  is.null(fit$null_delta),
  is.data.frame(fit$posterior_gene_estimates)
)

default_grid <- effect_size_grid(prevalence = 0.02)
wide_grid <- effect_size_grid(prevalence = 0.02, max_effect_size = 100)
stopifnot(
  length(default_grid) == 10,
  isTRUE(all.equal(exp(tail(default_grid, 1)), 50)),
  isTRUE(all.equal(exp(tail(wide_grid, 1)), 100))
)

reference <- burdenmle_gene_reference()
features <- loeuf_features(reference$gene_id[1:25], reference)
strata <- burdenmle_loeuf_strata()
stopifnot(
  nrow(reference) == 17395,
  nrow(features) == 25,
  ncol(features) == 6,
  all(rowSums(features) == 1),
  nrow(strata) == 6,
  !any(grepl("count|variant|sample", names(reference), ignore.case = TRUE))
)

bad <- dat[, setdiff(names(dat), "case_rate")]
input_error <- try(
  BurdenMLE_DN(
    bad,
    prevalence = 0.02,
    bootstrap = FALSE,
    null_sim = FALSE
  ),
  silent = TRUE
)
stopifnot(inherits(input_error, "try-error"))
