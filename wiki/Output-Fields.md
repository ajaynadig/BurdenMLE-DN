# Output fields

`BurdenMLE_DN()` returns a `BurdenMLEDN_fit` list. The most useful fields are:

| Field | Contents |
|---|---|
| `delta` | Fitted component weights; rows are annotation strata |
| `component_endpoints` | Upper log-rate-ratio endpoint for each component |
| `ll` | Maximized log likelihood |
| `optimizer_elapsed` | Elapsed time for the full fit and, when run, bootstrap |
| `mutvar_output$total_mutvar` | Total observed-scale mutational variance |
| `mutvar_output$mutvar_CI` | Bootstrap 95% interval for total mutational variance |
| `mutvar_output$annot_mutvar` | Annotation-stratum mutational variance |
| `mutvar_output$frac_mutvar` | Fraction of mutational variance by stratum |
| `mutvar_output$enrichment` | Mutational-variance enrichment by stratum |
| `penetrance$effective_penetrance` | Effective penetrance |
| `penetrance$effective_penetrance_CI` | Bootstrap 95% interval |
| `posterior_gene_estimates` | Optional gene-level naive and posterior summaries |
| `bootstrap_output` | Bootstrap weights and resampling indices |
| `null_delta` | Mixture weights from parametric null simulations |

Use `summary(fit)` for a compact collection and `mixture_weights(fit)` for the
fitted weights alone.
