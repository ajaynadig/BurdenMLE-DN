# Bundled reference files

`synthetic_example.csv` is entirely simulated.

`burdenmle_gene_reference.csv` contains only gene-level reference annotations:
Ensembl and gene identifiers, chromosome, LOEUF, the six baseline annotation
strata, predicted coding SNV mutation rates, and the gnomAD v4.1 synonymous
calibration factor used in the manuscript analysis. The `corrected_mu_*`
columns are the products of the predicted rates and that calibration factor.

`loeuf_strata_summary.csv` summarizes the six annotation strata. The first
LOEUF quintile is split at the global 80th percentile of corrected cumulative
coding SNV mutation rate; the remaining quintiles are not split.

These files contain no cohort observations, variant records, case or control
counts, sample identifiers, phenotypes, or fitted model estimates. PTV rates
in the table cover the indicated SNV class; users analyzing additional variant
types such as indels should supply mutation rates or a scaling procedure
appropriate to their own calling and annotation pipeline.
