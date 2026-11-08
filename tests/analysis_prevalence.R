# Run from the repository root (no governed data needed).
source("analysis/scripts/analysis_prevalence.R")
strata <- c("Male DDID", "Male Non-DDID", "Female DDID", "Female Non-DDID")
model_names <- c(paste(strata, "GeneDx PTV"), "Combined ASC PTV",
                 "Combined SPARK PTV", "Combined Probands PTV")
base <- c(.043 * .369, .043 * (1 - .369), .0114 * .421,
          .0114 * (1 - .421), rep(.0276, 3))
dat <- list(loop_vars = list(names = model_names, prevalences = base),
            prev_factors = c(.0221666598 / .0276, .01 / .0276, 1))
scenarios <- lapply(model_names, function(n) analysis_prevalence(dat, n))
stopifnot(identical(vapply(scenarios, `[[`, integer(1), "factor_index"),
                    c(2L, 2L, 2L, 2L, 3L, 3L, 1L)))
stopifnot(isTRUE(all.equal(vapply(scenarios, `[[`, numeric(1), "prevalence"),
                          c(base[1:4] * .01 / .0276, .0276, .0276, .0221666598))))
for (i in seq_along(scenarios)) {
  p <- scenarios[[i]]$prevalence
  assert_fit_prevalence(list(prevalence = p), p, model_names[i])
  stopifnot(inherits(try(assert_fit_prevalence(list(prevalence = p * 2.76),
                                              p, model_names[i]), silent = TRUE), "try-error"))
}
cat("Analysis prevalence tests passed.\n")
