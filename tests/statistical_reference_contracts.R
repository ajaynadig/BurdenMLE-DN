library(BurdenMLEDN)

expect_error <- function(code) {
  result <- try(force(code), silent = TRUE)
  stopifnot(inherits(result, "try-error"))
}

quiet_fit <- function(args) {
  fit <- NULL
  invisible(capture.output(fit <- do.call(BurdenMLE_DN, args)))
  fit
}

# Gene-level Poisson probabilities use the inclusive upper tail P(X >= x).
poisson_data <- data.frame(
  case_count = 0:2,
  expected_count = rep(0.7, 3),
  row.names = paste0("p", 0:2)
)
poisson_fit <- quiet_fit(list(
  input_data = poisson_data,
  component_endpoints = c(0, log(2)),
  bootstrap = FALSE,
  null_sim = FALSE,
  mutvar_est = FALSE,
  estimate_effective_penetrance = FALSE
))
expected_poisson <- c(
  1,
  1 - exp(-0.7),
  1 - exp(-0.7) * (1 + 0.7)
)
observed_poisson <- poisson_fit$posterior_gene_estimates$Estimate_Poisson_P
stopifnot(isTRUE(all.equal(observed_poisson, expected_poisson, tolerance = 1e-14)))

# Finite-null probabilities include ties, never equal zero, and are no smaller
# than the former strict empirical tail.
finite_null_p <- BurdenMLEDN:::inclusive_finite_null_pvalue
null_cases <- list(
  below = list(values = c(0, 1, 2), observed = 3, expected = 1 / 4),
  tied = list(values = c(1, 2, 2), observed = 2, expected = 3 / 4),
  above = list(values = c(3, 4, 5), observed = 2, expected = 1)
)
for (case in null_cases) {
  corrected <- finite_null_p(case$values, case$observed)
  historical <- mean(case$values > case$observed)
  stopifnot(
    identical(corrected, case$expected),
    corrected > 0,
    corrected >= historical
  )
}

# The forecasting script's integer-count expression is the same inclusive tail.
counts <- 0:5
lambda <- 0.8
stopifnot(isTRUE(all.equal(
  ppois(counts - 1, lambda, lower.tail = FALSE),
  ppois(counts - 0.001, lambda, lower.tail = FALSE),
  tolerance = 0
)))

# A two-gene fixture makes each gene's component and midpoint effect exact.
make_penetrance_fixture <- function(case_rate) {
  input <- data.frame(
    case_count = c(0, 0),
    case_rate = case_rate,
    N = c(100, 100),
    row.names = c("low", "high")
  )
  genetic_data <- BurdenMLEDN:::process_data_trio(input)
  features <- diag(2)
  dimnames(features) <- list(rownames(input), c("low", "high"))
  fit <- BurdenMLEDN:::initialize_model(
    BurdenMLEDN:::poisson_uniform_likelihood,
    genetic_data,
    component_endpoints = c(log(2), log(8)),
    features = features,
    grid_size = 1L
  )
  fit$delta <- diag(2)
  dimnames(fit$delta) <- list(colnames(features), NULL)
  fit$posterior_gene_estimates <- data.frame(
    Case_Count = genetic_data$case_count,
    Expected_Count = genetic_data$expected_count,
    row.names = rownames(genetic_data)
  )
  fit$prevalence <- 0.02
  class(fit) <- c("BurdenMLEDN_fit", "list")
  list(fit = fit, input = input, genetic_data = genetic_data)
}

heterogeneous <- make_penetrance_fixture(c(0.001, 0.004))
gamma <- sqrt(c(2, 8))
numerator <- (gamma - 1) * gamma
denominator <- gamma - 1
gene_average <- 0.02 * mean(numerator) / mean(denominator)
mutation_weighted <- 0.02 *
  sum(heterogeneous$input$case_rate * numerator) /
  sum(heterogeneous$input$case_rate * denominator)
internal_weighted <- BurdenMLEDN:::effective_penetrance_func(
  heterogeneous$fit, heterogeneous$genetic_data, 0.02
)
stopifnot(
  isTRUE(all.equal(internal_weighted, mutation_weighted, tolerance = 1e-14)),
  internal_weighted > gene_average
)

symmetric <- make_penetrance_fixture(c(0.002, 0.002))
stopifnot(isTRUE(all.equal(
  BurdenMLEDN:::effective_penetrance_func(
    symmetric$fit, symmetric$genetic_data, 0.02
  ),
  gene_average,
  tolerance = 1e-14
)))

expected_only <- heterogeneous$input[c("case_count")]
expected_only$expected_count <- heterogeneous$genetic_data$expected_count
expect_error(BurdenMLEDN:::effective_penetrance_func(
  heterogeneous$fit, expected_only, 0.02
))

# The default fit exposes the mutation-weighted estimand and therefore requires
# explicit mutation rates.
default_penetrance_data <- data.frame(
  case_count = c(0, 1, 0, 2),
  case_rate = c(0.001, 0.002, 0.0015, 0.0025),
  N = rep(100, 4),
  row.names = paste0("d", 1:4)
)
expected_count_only <- data.frame(
  case_count = default_penetrance_data$case_count,
  expected_count = 2 * default_penetrance_data$N *
    default_penetrance_data$case_rate,
  row.names = rownames(default_penetrance_data)
)
expect_error(quiet_fit(list(
  input_data = expected_count_only,
  component_endpoints = c(0, log(2), log(5)),
  prevalence = 0.02,
  bootstrap = FALSE,
  null_sim = FALSE,
  mutvar_est = FALSE,
  estimate_effective_penetrance = TRUE
)))
default_penetrance_fit <- quiet_fit(list(
  input_data = default_penetrance_data,
  component_endpoints = c(0, log(2), log(5)),
  prevalence = 0.02,
  bootstrap = FALSE,
  null_sim = FALSE,
  mutvar_est = FALSE,
  estimate_effective_penetrance = TRUE
))
stopifnot(identical(
  names(default_penetrance_fit$penetrance),
  c("effective_penetrance", "effective_penetrance_CI")
))
default_genetic_data <- BurdenMLEDN:::process_data_trio(default_penetrance_data)
stopifnot(isTRUE(all.equal(
  default_penetrance_fit$penetrance$effective_penetrance,
  BurdenMLEDN:::effective_penetrance_func(
    default_penetrance_fit, default_genetic_data, 0.02
  ),
  tolerance = 1e-14
)))

# The stored bootstrap interval uses the resampled genes' corresponding
# mutation weights in every replicate.
penetrance_bootstrap_samples <- matrix(
  c(
    "d1", "d2", "d3", "d4",
    "d1", "d1", "d2", "d4",
    "d2", "d3", "d3", "d4",
    "d1", "d2", "d4", "d4"
  ),
  nrow = 4,
  dimnames = list(NULL, paste0("bootstrap_", 1:4))
)
bootstrap_penetrance_fit <- quiet_fit(list(
  input_data = default_penetrance_data,
  component_endpoints = c(0, log(2), log(5)),
  prevalence = 0.02,
  bootstrap = TRUE,
  bootstrap_samples = penetrance_bootstrap_samples,
  n_boot = 4,
  null_sim = FALSE,
  mutvar_est = FALSE,
  estimate_effective_penetrance = TRUE
))
bootstrap_genetic_data <- BurdenMLEDN:::process_data_trio(
  default_penetrance_data
)
bootstrap_penetrance_estimates <- vapply(seq_len(4), function(iteration) {
  replicate <- BurdenMLEDN:::reconstruct_bootstrap_replicate(
    bootstrap_penetrance_fit,
    bootstrap_genetic_data,
    iteration
  )
  BurdenMLEDN:::effective_penetrance_func(
    replicate$model,
    replicate$genetic_data,
    prevalence = 0.02
  )
}, numeric(1))
stopifnot(isTRUE(all.equal(
  unname(bootstrap_penetrance_fit$penetrance$effective_penetrance_CI),
  unname(quantile(bootstrap_penetrance_estimates, c(0.025, 0.975))),
  tolerance = 1e-14
)))

# The bundled resources satisfy logical invariants, and targeted corruptions do
# not pass the same private validators. No file fingerprint is asserted.
reference <- burdenmle_gene_reference()
strata <- burdenmle_loeuf_strata()
validate_reference <- BurdenMLEDN:::validate_burdenmle_gene_reference
validate_strata <- BurdenMLEDN:::validate_burdenmle_loeuf_strata
stopifnot(
  is.data.frame(validate_reference(reference)),
  is.data.frame(validate_strata(strata, reference)),
  identical(rownames(loeuf_features(reference$gene_id[1:10])),
            reference$gene_id[1:10]),
  all(rowSums(loeuf_features(reference$gene_id[1:10])) == 1)
)

bad_reference <- reference[-1]
expect_error(validate_reference(bad_reference))
bad_reference <- reference
bad_reference$gene_id[2] <- bad_reference$gene_id[1]
expect_error(validate_reference(bad_reference))
bad_reference <- reference
bad_reference$mu_snp_PTV[1] <- -1
expect_error(validate_reference(bad_reference))
bad_reference <- reference
bad_reference$corrected_mu_snp_PTV[1] <- 0
expect_error(validate_reference(bad_reference))
bad_reference <- reference[c(2, 1, seq.int(3, nrow(reference))), ]
expect_error(validate_reference(bad_reference))
bad_reference <- reference
bad_reference$analysis_stratum[1] <- "unknown"
expect_error(validate_reference(bad_reference))
bad_reference <- reference
interior_index <- function(stratum) {
  rows <- which(reference$analysis_stratum == stratum)
  rows[reference$loeuf[rows] > min(reference$loeuf[rows]) &
    reference$loeuf[rows] < max(reference$loeuf[rows])][1]
}
swap <- c(interior_index("LOEUF1_mu1"), interior_index("LOEUF1_mu2"))
bad_reference$analysis_stratum[swap] <-
  rev(bad_reference$analysis_stratum[swap])
expect_error(validate_reference(bad_reference))

bad_strata <- strata
bad_strata$genes[1] <- bad_strata$genes[1] + 1L
expect_error(validate_strata(bad_strata, reference))
bad_strata <- strata[rev(seq_len(nrow(strata))), ]
expect_error(validate_strata(bad_strata, reference))
bad_strata <- strata
bad_strata$loeuf_max[1] <- bad_strata$loeuf_max[1] + 0.1
expect_error(validate_strata(bad_strata, reference))

# The minimal custom-reference seam of loeuf_features() remains valid.
custom_reference <- data.frame(
  gene_id = c("custom1", "custom2"),
  analysis_stratum = c("LOEUF2", "LOEUF5")
)
custom_features <- loeuf_features(custom_reference$gene_id, custom_reference)
stopifnot(
  identical(rownames(custom_features), custom_reference$gene_id),
  all(rowSums(custom_features) == 1)
)
