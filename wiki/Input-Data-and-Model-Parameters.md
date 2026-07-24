# Input data and model parameters

## Input table

`BurdenMLE_DN()` accepts a data frame with one row per gene, unique gene row
names, and three required columns:

| Field | Meaning |
|---|---|
| `case_count` | Nonnegative integer number of observed de novo variants |
| `case_rate` | Per-haploid gene-level mutation rate for the analyzed class |
| `N` | Number of affected offspring; repeat the same value for every gene |

The expected count under no association is computed as
`2 * N * case_rate`. Rows with missing or non-finite values are rejected.

For stratified models, `features` is a one-hot numeric matrix with one row per
gene. Its row names and order must match the input table. Each column defines
an annotation stratum with its own mixture weights.

## Principal model parameters

| Argument | Purpose | Default |
|---|---|---|
| `prevalence` | Population prevalence used for effect grid and derived estimates | required |
| `no_cpts` | Number of mixture components | 10 |
| `component_endpoints` | Optional explicit upper endpoints on the log-rate-ratio scale | prevalence-dependent grid |
| `grid_size` | Integration points per uniform component | 10 |
| `bootstrap`, `n_boot` | Gene-level nonparametric bootstrap and replicate count | `TRUE`, 100 |
| `null_sim`, `n_null` | Parametric null simulation and replicate count | `TRUE`, 100 |
| `estimate_posteriors` | Gene-level posterior summaries | `FALSE` |
| `optimizer` | Mixture optimizer | `"mixsqp"` |
| `mixsqp_control` | Overrides passed to `mixsqp::mixsqp()` | package defaults |

The MixSQP defaults use no preliminary iterations, retain numerically small
components (`zero.threshold.solution = 0`), and set the active-set iteration
allowance to `max(20, 2 * number_of_components)`.
