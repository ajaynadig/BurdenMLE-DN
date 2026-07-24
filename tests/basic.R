library(BurdenMLEDN)

path <- system.file("extdata", "synthetic_example.csv", package = "BurdenMLEDN")
dat <- read.csv(path, row.names = "gene")

set.seed(1)
fit <- BurdenMLE_DN(
  dat,
  prevalence = 0.02,
  no_cpts = 6,
  bootstrap = FALSE,
  null_sim = FALSE,
  mutvar_est = FALSE,
  estimate_effective_penetrance = FALSE
)

stopifnot(
  inherits(fit, "BurdenMLEDN_fit"),
  identical(fit$optimizer, "mixsqp"),
  isTRUE(all.equal(sum(fit$delta), 1, tolerance = 1e-6)),
  is.finite(fit$ll)
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
