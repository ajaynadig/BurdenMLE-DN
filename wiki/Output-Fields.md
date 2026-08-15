# Output fields

`BurdenMLE_DN()` returns a `BurdenMLEDN_fit` list. All stratum-indexed outputs
retain the feature-column names and order. The most useful fields are:

| Field | Availability | Description |
|---|---|---|
| `delta` | Always | Fitted mixture probabilities. Rows correspond to annotation strata and columns to effect-size components. |
| `component_endpoints` | Always | Upper log-rate-ratio endpoint of each uniform mixture component. Exponentiate these values to recover rate-ratio endpoints. |
| `ll` | When `return_likelihood = TRUE` | Maximized log likelihood, useful for fit diagnostics and optimizer comparisons on identical data and grids. |
| `fit_status` | Always | Backend-neutral full-fit record containing the final log likelihood and weights, convergence/usability flags, stable status code, backend message, and iteration count. |
| `uncertainty_reliable` | Always | `TRUE` when every requested bootstrap/null fit converged, `FALSE` when a finite usable nonconverged replicate contributed, and `NA` when no uncertainty stage was requested. |
| `optimizer_elapsed` | Always | Elapsed seconds for the full fit and, when requested, the complete bootstrap fitting stage. |
| `input_summary` | Always | Input gene/count totals, sample size (`NA` when `N` was not supplied), and whether `case_rate` is available. |
| `mutvar_output$total_mutvar` | When `mutvar_est = TRUE` | Total observed-scale coding mutational variance summed across annotation strata using deterministic analytic moments of the continuous-uniform components. |
| `mutvar_output$mutvar_CI` | With bootstrap and mutational variance | Percentile-bootstrap 95% interval for total mutational variance. |
| `mutvar_output$annot_mutvar` | When `mutvar_est = TRUE` | Absolute contribution to mutational variance from each annotation stratum. |
| `mutvar_output$annot_mutvar_CI` | With bootstrap and mutational variance | Stratum-specific percentile-bootstrap 95% intervals for absolute mutational variance. |
| `mutvar_output$frac_mutvar` | When total mutational variance is positive | Fraction of total mutational variance assigned to each annotation stratum. |
| `mutvar_output$fracmutvar_CI` | With a defined bootstrap fraction | Percentile-bootstrap 95% intervals for the stratum fractions. |
| `mutvar_output$frac_expected` | When `mutvar_est = TRUE` | Fraction of the total reference mutation rate contributed by each annotation stratum. |
| `mutvar_output$enrichment` | When total mutational variance is positive | `frac_mutvar / frac_expected`; values above one indicate more mutational variance than expected from mutation rate alone. |
| `mutvar_output$enrich_CI` | With a defined bootstrap enrichment | Percentile-bootstrap 95% intervals for enrichment. |
| `penetrance$effective_penetrance` | When `estimate_effective_penetrance = TRUE` | Gene-average penetrance among the excess risk attributable to the modeled variant class. This remains the package's default effective-penetrance estimand. |
| `penetrance$effective_penetrance_CI` | With bootstrap and penetrance | Percentile-bootstrap 95% interval for gene-average effective penetrance. |
| `posterior_gene_estimates` | When `estimate_posteriors = TRUE` | Observed/expected rate ratios, inclusive one-sided Poisson probabilities `P(X >= observed)`, and posterior-mean gene rate ratios. |
| `bootstrap_output` | When `bootstrap = TRUE` | Named replicate fit records and compatibility weight matrices, canonical sampled gene IDs, resolved local indices, and replay seeds. |
| `null_output` | Only when `null_sim = TRUE` | Named null-replicate fit records, compatibility weight matrices, and replay seeds. |
| `null_delta` | Only when `null_sim = TRUE` | Compatibility view of the null mixture weights; not needed for point estimates or bootstrap confidence intervals. |

Use `summary(fit)` for a compact collection and `mixture_weights(fit)` for the
fitted weights alone. Full-data nonconvergence is an error. A finite,
simplex-valid nonconverged uncertainty replicate remains in the interval and
sets `uncertainty_reliable = FALSE`; a nonfinite or invalid replicate stops the
uncertainty stage.

Call `mutation_weighted_effective_penetrance(fit, input_data)` explicitly to
calculate the distinct mutation-weighted estimand. It is not computed or used
by default.

When null simulation is requested, `mutvar_output$mutvar_p` includes null
values tied with or above the observed mutational variance and uses the
finite-null plus-one correction. Its minimum possible value is
`1 / (n_null + 1)`.
