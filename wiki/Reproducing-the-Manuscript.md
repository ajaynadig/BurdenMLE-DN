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
4. Copy `analysis/workflow/manuscript.json`, set its absolute `root` and run identifier, and inspect the plan before executing:

```bash
bash analysis/reproduce_study.sh --config /path/to/run.json
bash analysis/reproduce_study.sh --config /path/to/run.json --backend local --run
```

The runner writes models, derived results, figures, tables, and logs beneath
the supplied analysis root. Local and SLURM execution share the same dependency
graph. Select targets with `--stages`; explicitly validate/reuse completed stages
with `--reuse`. The full study, cluster configuration, and reporting-only imports are documented in
[`analysis/README.md`](../analysis/README.md).

Every model-producing mode also writes
`outputs/data/model_manifest_<mode>_<run-date>.rds`. Downstream stages load the
exact artifact named in that manifest and validate its saved objects, gene
order, full-fit status, and uncertainty status; they never choose a model by
sorting filenames. Reporting-only rebuilds explicitly import validated approved
results into a new root; they do not imply fresh model fitting. Direct R-script
calls still accept `--model-manifest`. Unreliable bootstrap uncertainty warns without
blocking a usable fit. Historical unmanifested files require an explicit
`--legacy-...-model-file` option and always warn because their modern status
and gene-universe contracts cannot be reconstructed.

This separation is deliberate: users can install and test the method with
synthetic data, while controlled study data remain outside Git and under the
governance of the flagship manuscript authors.
