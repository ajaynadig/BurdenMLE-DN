# Manuscript reproduction

This directory contains the scripts used for the BurdenMLE-DN manuscript.
Model fitting, summaries, figures, tables, and the manuscript-estimate
checklist are selected through `reproduce_study.sh`.

The manuscript cohort inputs are **not distributed in this repository**.
They include embargoed gene-count and variant tables governed by the
originating study authors. Authorized analysts should create an analysis root
with the directory structure in `input_manifest.tsv`, then run:

```bash
BURDENMLEDN_ANALYSIS_ROOT=/absolute/path/to/authorized_analysis_root \
RUN_ALL=true \
bash analysis/reproduce_study.sh
```

If `BURDENMLEDN_ANALYSIS_ROOT` is omitted, the runner uses this `analysis/`
directory. Its `inputs/`, `outputs/`, and `logs/` directories are ignored by
Git.

All stages are opt-in unless `RUN_ALL=true`. Individual stages can be selected
with the `RUN_*` variables declared near the top of `reproduce_study.sh`.
`RUN_DATE`, seeds, bootstrap counts, and R executable can also be set as
environment variables.

The scripts are released for transparency and exact reproduction. The
installable package and synthetic example do not require manuscript inputs.
