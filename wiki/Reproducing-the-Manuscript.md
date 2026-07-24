# Reproducing the manuscript

The complete analysis scripts are in [`analysis/`](../analysis). The repository
does not distribute the manuscript cohorts' variant, count, pedigree, or
phenotype tables.

Authorized analysts should:

1. Install this package and the analysis dependencies.
2. Create an analysis root following
   [`analysis/input_manifest.tsv`](../analysis/input_manifest.tsv).
3. Obtain each restricted input from its originating study authors and each
   public reference from the cited source.
4. Run:

```bash
BURDENMLEDN_ANALYSIS_ROOT=/absolute/path/to/authorized_analysis_root \
RUN_ALL=true \
bash analysis/reproduce_study.sh
```

The runner writes models, derived results, figures, tables, and logs beneath
the supplied analysis root. Every stage is separately selectable through the
`RUN_*` variables documented in
[`analysis/README.md`](../analysis/README.md).

This separation is deliberate: users can install and test the method with
synthetic data, while controlled study data remain outside Git and under the
governance of the flagship manuscript authors.
