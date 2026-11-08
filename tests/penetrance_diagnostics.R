library(BurdenMLEDN)

capture_warnings <- function(code) {
  messages <- character()
  value <- withCallingHandlers(force(code), warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, messages = messages)
}

# Diagnostics do not change finite values or the percentile interval policy.
diagnostic <- BurdenMLEDN:::warn_effective_penetrance
interval <- BurdenMLEDN:::effective_penetrance_interval
stopifnot(length(capture_warnings(diagnostic(c(0, 0.5, 1)))$messages) == 0L)
for (values in list(c(0.1, 0.4, 0.8), c(-0.1, 0.5, 1.5))) {
  result <- capture_warnings(interval(values))
  stopifnot(identical(result$value, quantile(values, c(0.025, 0.975))))
  stopifnot(length(result$messages) == as.integer(any(values < 0 | values > 1)))
}
invalid <- capture_warnings(interval(c(0.5, 1.2, NA_real_, Inf)))
stopifnot(
  identical(invalid$value, c("2.5%" = NA_real_, "97.5%" = NA_real_)),
  length(invalid$messages) == 1L,
  grepl("1 value(s) outside [0,1] and 2 non-finite", invalid$messages, fixed = TRUE),
  grepl("CI will not be calculated", invalid$messages, fixed = TRUE)
)

input <- data.frame(case_count = c(0, 0), case_rate = c(0.001, 0.004),
                    N = c(100, 100), row.names = c("a", "b"))
data <- BurdenMLEDN:::process_data_trio(input)
features <- matrix(1, 2, 1, dimnames = list(rownames(input), "all"))
fit <- BurdenMLEDN:::initialize_model(
  BurdenMLEDN:::poisson_uniform_likelihood, data,
  component_endpoints = c(0, log(4)), features = features, grid_size = 1L
)
fit$delta <- matrix(c(0, 1), 1, 2, dimnames = list("all", NULL))
point <- BurdenMLEDN:::effective_penetrance_func
moments <- BurdenMLEDN:::precompute_effective_penetrance_moments(fit, data)
ordinary <- capture_warnings(point(fit, data, 0.02, moments = moments))
large <- capture_warnings(point(fit, data, 0.8, moments = moments))
stopifnot(
  isTRUE(all.equal(ordinary$value, 0.04)), length(ordinary$messages) == 0L,
  isTRUE(all.equal(large$value, 1.6)), length(large$messages) == 1L
)
# Non-finite numerator and exact-null denominator retain their existing results.
bad_moments <- moments
bad_moments$numerator[,] <- Inf
nonfinite <- capture_warnings(point(fit, data, 0.02, moments = bad_moments))
stopifnot(!is.finite(nonfinite$value), length(nonfinite$messages) == 1L,
          identical(nonfinite$value, BurdenMLEDN:::effective_penetrance_from_moments(
            fit, bad_moments, 0.02, input$case_rate)))
fit$delta[,] <- c(1, 0)
null <- capture_warnings(point(fit, data, 0.02, moments = moments))
stopifnot(is.na(null$value), length(null$messages) == 1L,
          grepl("near the null", null$messages, fixed = TRUE))
# No new near-null threshold: tiny non-null weight still evaluates the ratio.
fit$delta[,] <- c(1 - 1e-12, 1e-12)
near_null <- capture_warnings(point(fit, data, 0.02, moments = moments))
stopifnot(isTRUE(all.equal(near_null$value, 0.04)),
          length(near_null$messages) == 0L)

# Full estimator emits one aggregate bootstrap warning, not one per replicate.
large_input <- input
large_input$case_count <- c(5, 8)
large_fit <- NULL
large_run <- capture_warnings(invisible(capture.output(
  large_fit <- BurdenMLE_DN(large_input, component_endpoints = c(0, log(4)),
                            prevalence = 0.8, n_boot = 4,
                            bootstrap_seed = 19, mutvar_est = FALSE)
)))
stopifnot(large_fit$penetrance$effective_penetrance > 1,
          all(large_fit$penetrance$effective_penetrance_CI > 1),
          sum(grepl("Bootstrap effective penetrance", large_run$messages,
                    fixed = TRUE)) == 1L)

reference_messages <- character()
reference <- withCallingHandlers(burdenmle_gene_reference(), message = function(m) {
  reference_messages <<- c(reference_messages, conditionMessage(m))
  invokeRestart("muffleMessage")
})
stopifnot(nrow(reference) == 17395L,
          any(grepl("Satterstrom", reference_messages, fixed = TRUE)),
          any(grepl("callable exposure", reference_messages, fixed = TRUE)))
