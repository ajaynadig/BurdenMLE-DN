# Manuscript reproduction

`reproduce_study.sh` runs one dependency graph locally or on SLURM. The package
and [public Kaplanis example](../inst/examples/kaplanis/README.md) do not require
the authorized manuscript inputs. The example is not the manuscript DD analysis.

## Prepare an environment and inputs

Use Python 3.9+ and R with the package dependencies plus the analysis dependencies
checked by `workflow/preflight.R` (including Bioconductor SummarizedExperiment,
ggplot2, patchwork, data.table, readr, dplyr, tidyr, scales, MetBrewer, openxlsx,
and digest). Preflight installs this checkout of BurdenMLEDN into the run's own
`R-library` and records installed package versions. It does not install missing
analysis dependencies automatically. For SLURM, provide the same environment on
every node, optionally through `slurm.image` and explicit `slurm.binds`.

Place authorized inputs in `<root>/inputs/` using `input_manifest.tsv`. Counts and
variant tables are governed by their originating studies and are not on GitHub.
The public CDC/NCHS `Nat2024us.zip` archive must also be extracted under
`inputs/demographics/` before launch. It uses Deflate64: use a compatible unzip
(for example macOS `/usr/bin/unzip`). Preflight verifies the extracted member's
size and hashes the inputs; compute jobs read the extracted text directly.

Copy `workflow/manuscript.json` outside the repository and set `root` to a new
absolute directory and `run_id` to a stable identifier. Set `bootstrap_draws` to
the approved directory of saved gene-ID draws to reproduce manuscript resamples;
leave it empty to generate new deterministic draws. Do not expect a seed alone
to reconstruct draws made by a different historical execution order. Saved draws
are validated against each model's ordered gene universe and bootstrap count.
They may contain restricted identifiers and are not distributed here.

The `manuscript-r2` profile fixes 100 gene bootstraps, 500 main simulation
replicates, 100 matched optimizer-comparison replicates, 25 replicates per batch,
and upper RR bounds 100/250. `analysis_prevalence.R` selects pooled, clinical,
or population scenarios from the original setup definitions: combined-sex
ASC/SPARK use 2.76%, GeneDx 1%, and pooled ASD remains unchanged. Stratum-specific
prevalences remain stratum-specific. Reporting verifies agreement with saved fits.

The CES-excluded upper-bound sweep fits ASD and DD PTV and Mis2 models at ten
bounds (40 tasks, each with 100 bootstraps). Within each disorder, the classes
share ordered genes and bootstrap draws. Likelihood, variance, and case fraction
combine PTV + Mis2; effective penetrance is PTV-only for both disorders.

## Run

```bash
# Print the entire graph without launching anything.
bash analysis/reproduce_study.sh --config /path/to/run.json

# Same stages and random inputs, executed sequentially on this machine:
bash analysis/reproduce_study.sh --config /path/to/run.json --backend local --run

# Run this command on a SLURM submission host (sbatch must be on PATH):
bash analysis/reproduce_study.sh --config /path/to/run.json --backend slurm --run
```

SLURM uses after-success dependencies and independent model, simulation and
sensitivity arrays. `slurm.array_limit` is null by default (no application-level
cap); scheduler/account limits still apply. `forecast_cores` controls deterministic
within-job forecasting parallelism. SLURM settings contain no lab-specific paths.
Local execution blocks until completion; SLURM execution records job IDs in
`workflow/submitted.json` and returns after submission, not after completion.

`--stages figure2,figure3` selects targets and their dependencies. Nothing is
implicitly reused. `--reuse models` verifies a completed models stage and stops
dependency traversal there. To rerun reporting from existing results, select the
required targets and explicitly reuse their completed dependencies. Existing
completed stages are not overwritten. Use a new root for changed code/settings.
Legacy RUN_ALL/RUN_* selectors are rejected, including flags set to false; unset
them and migrate to explicit `--stages`/`--reuse`. This avoids silently broadening
old partial-run commands. RUN_DATE remains a supported run identifier when using
the default configuration. No environment selector implicitly authorizes launch.

## Rebuild reporting from approved fits

This is a different operation from a fresh full study. Import into a new root
with the same run_id as the approved source:

```bash
python3 analysis/workflow/import_results.py --config /path/to/report.json \
  --source /path/to/approved/run
bash analysis/reproduce_study.sh --config /path/to/report.json --run \
  --reuse models,simulations,upper,korean,public_example
```

Provide the authorized inputs and an installed analysis environment in the new
root as above. Import runs preflight to record raw-input, draw and environment
identities, then copies selected fitted/result artifacts, records source
hashes, relocates model manifests to the copies, and validates their schema and
fit records. It does not copy empirical reporting summaries or forecasts. The
import manifest explicitly labels this as a reporting rebuild, not fresh fitting.

## Outputs, checks, and restart

`outputs/` holds analysis products; `submission/` holds numbered figures, ten
supplementary TSVs and the combined workbook, numbering maps, and the 214-row
manuscript checklist. `provenance/` records input/environment/output identities;
`workflow/` records configuration, code hashes, stage completion and submitted
jobs. Logs are per stage/task. Historical R2 comparisons are optional external
audit work, not a dependency of ordinary reproduction.

Completion markers are written only after success. Reuse checks recorded output
hashes and configuration/code identity. Failed workers do not release dependent
SLURM jobs. Inspect failures and any stale locks before retrying; do not remove a
lock while its worker is alive. A submission ledger prevents accidental duplicate
submission. Archive it only after checking all previously submitted jobs have
finished or been cancelled. A stage that failed after producing partial outputs
may require a fresh root; the runner favors refusing ambiguous overwrites.

The generated workbook is read back and compared to its TSVs. Automated checks
do not replace visual inspection or manuscript cross-reference review. Freeze
only after those checks and pin the code commit plus configuration, input/draw
hashes, environment, and output manifest. Inputs must remain available under
their access rules; external symlinks are not an independent archive.

Developer tests: run `python3 tests/workflow.py` and the R scripts in `tests/`
from a source checkout with the package installed. Native SLURM submission and
the authorized-data reporting rebuild require their respective environments.
