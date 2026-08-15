# Input data and model parameters

## Input table

`BurdenMLE_DN()` accepts a data frame with at least two genes, explicit,
unique, nonempty, non-missing gene row names, and three required columns:

| Field | Meaning |
|---|---|
| `case_count` | Nonnegative integer number of observed de novo variants |
| `case_rate` | Per-haploid gene-level mutation rate for the analyzed class |
| `N` | Number of affected offspring; repeat the same value for every gene |

The expected count under no association is computed as
`2 * N * case_rate`. Rows with missing or non-finite values are rejected.

For stratified models, `features` must be a finite one-hot numeric matrix with
one row per gene. Its unique, nonempty, non-missing row names must be identical
to the input gene names in the same order; reordered rows fail rather than
being silently realigned. Columns must have unique, nonempty, non-missing
names, every row must assign its gene to exactly one stratum, and every stratum
must contain a gene. The same contract applies to MixSQP and EM. Omitting
`features` fits one stratum named `all_genes`.

## Principal model parameters

| Argument | Purpose | Default |
|---|---|---|
| `prevalence` | Population prevalence used for effect grid and derived estimates | required |
| `no_cpts` | Number of mixture components | 10 |
| `component_endpoints` | Optional finite, unique, strictly increasing upper endpoints on the log-rate-ratio scale; takes precedence over `no_cpts` | prevalence-dependent grid |
| `grid_size` | Integration points per uniform component | 10 |
| `bootstrap`, `n_boot` | Gene-level nonparametric bootstrap and replicate count | `TRUE`, 100 |
| `null_sim`, `n_null` | Optional parametric null simulation and replicate count | `FALSE`, 100 when enabled |
| `estimate_posteriors` | Gene-level posterior summaries | `TRUE` |
| `optimizer` | Mixture optimizer | `"mixsqp"` |
| `mixsqp_control` | Overrides passed to `mixsqp::mixsqp()` | package defaults |

The MixSQP defaults use no preliminary iterations, retain numerically small
components (`zero.threshold.solution = 0`), and set the active-set iteration
allowance to `max(20, 2 * number_of_components)`.

## Choosing the effect-size grid

The default grid runs from rate ratio 1 to `1 / prevalence`. This is a
reasonable automatic choice when prevalence is well specified. If prevalence
may be misspecified, or very large effects are plausible, the default can
truncate meaningful support. In that setting, specify a conservatively large
upper bound:

```r
component_endpoints <- effect_size_grid(
  prevalence = 0.02,
  max_effect_size = 100
)
```

The manuscript used 10 components with maximum rate ratios of 100 for autism
and 250 for developmental disorders. These are analysis choices, not
disease-specific package defaults.

## Bundled gene reference

`burdenmle_gene_reference()` returns 17,395 autosomal genes with complete
reference information. It contains no cohort observations or variant records.
The fields include:

- Ensembl and gene identifiers, chromosome, and LOEUF;
- the five LOEUF quintiles used in the main analysis;
- a six-level `analysis_stratum`, splitting the most constrained LOEUF
  quintile by corrected cumulative coding SNV mutation rate;
- predicted gene-level PTV, Mis2, Mis1, Mis0, and synonymous SNV mutation
  rates;
- the gnomAD v4.1 synonymous calibration factor and corrected mutation rates.

The corresponding feature matrix is constructed directly:

```r
reference <- burdenmle_gene_reference()
burdenmle_loeuf_strata()
features <- loeuf_features(reference$gene_id)
```

For a new cohort, match its observed gene counts to `gene_id`, choose the
appropriate `corrected_mu_snp_*` column as `case_rate`, and use the same gene
order for `loeuf_features()`. PTV rates cover the indicated SNV class; analyses
that also include indels must supply a suitable indel mutation rate or scaling
procedure for their own calling pipeline.

## Manuscript settings versus package defaults

| Setting | Package default | Manuscript analysis |
|---|---|---|
| Effect grid | Maximum RR `1 / prevalence` | Maximum RR 100 for autism; 250 for developmental disorders |
| Features | One mixture if omitted | Six LOEUF/mutation-rate strata |
| Bootstrap | 100 gene-resampling replicates | 100 replicates |
| Gene posteriors | Calculated | Calculated |
| Null simulation | Off | Off |
| Optimizer | MixSQP | MixSQP |
| Components / integration grid | 10 / 10 | 10 / 10 |

The main models also reused identical bootstrap resampling indices across
related fits to make comparisons less sensitive to Monte Carlo variation.
