# Output fields

`BurdenMLE_DN()` returns a `BurdenMLEDN_fit` list. All stratum-indexed outputs
retain the feature-column names and order. The most useful fields are:

| Field | Availability | Description |
|---|---|---|
| `delta` | Always | Fitted mixture probabilities. Rows correspond to annotation strata and columns to effect-size components. |
| `component_endpoints` | Always | Upper log-rate-ratio endpoint of each uniform mixture component. Exponentiate these values to recover rate-ratio endpoints. |
| `ll` | When `return_likelihood = TRUE` | Maximized log likelihood, useful for fit diagnostics and optimizer comparisons on identical data and grids. |
| `optimizer_elapsed` | Always | Elapsed seconds for the full fit and, when requested, the complete bootstrap fitting stage. |
| `input_summary` | Always | Input gene/count totals, sample size (`NA` when `N` was not supplied), and whether `case_rate` is available. |
| `mutvar_output$total_mutvar` | When `mutvar_est = TRUE` | Total observed-scale coding mutational variance summed across annotation strata. |
| `mutvar_output$mutvar_CI` | With bootstrap and mutational variance | Percentile-bootstrap 95% interval for total mutational variance. |
| `mutvar_output$annot_mutvar` | When `mutvar_est = TRUE` | Absolute contribution to mutational variance from each annotation stratum. |
| `mutvar_output$annot_mutvar_CI` | With bootstrap and mutational variance | Stratum-specific percentile-bootstrap 95% intervals for absolute mutational variance. |
| `mutvar_output$frac_mutvar` | When total mutational variance is positive | Fraction of total mutational variance assigned to each annotation stratum. |
| `mutvar_output$fracmutvar_CI` | With a defined bootstrap fraction | Percentile-bootstrap 95% intervals for the stratum fractions. |
| `mutvar_output$frac_expected` | When `mutvar_est = TRUE` | Fraction of the total reference mutation rate contributed by each annotation stratum. |
| `mutvar_output$enrichment` | When total mutational variance is positive | `frac_mutvar / frac_expected`; values above one indicate more mutational variance than expected from mutation rate alone. |
| `mutvar_output$enrich_CI` | With a defined bootstrap enrichment | Percentile-bootstrap 95% intervals for enrichment. |
| `penetrance$effective_penetrance` | When `estimate_effective_penetrance = TRUE` | Mutation-weighted average penetrance among the excess risk attributable to the modeled variant class. |
| `penetrance$effective_penetrance_CI` | With bootstrap and penetrance | Percentile-bootstrap 95% interval for effective penetrance. |
| `posterior_gene_estimates` | When `estimate_posteriors = TRUE` | Observed/expected rate ratios, Poisson tail probabilities, and posterior-mean gene rate ratios. |
| `bootstrap_output` | When `bootstrap = TRUE` | Fitted bootstrap mixture weights and the gene-resampling indices used in each replicate. |
| `null_delta` | Only when `null_sim = TRUE` | Mixture weights from optional parametric null datasets; not needed for point estimates or bootstrap confidence intervals. |

Use `summary(fit)` for a compact collection and `mixture_weights(fit)` for the
fitted weights alone.
