# Public Kaplanis developmental-disorder example

This short interactive script downloads Supplementary Table 1 from Kaplanis
et al. (*Nature* 2020; doi:10.1038/s41586-020-2832-5), counts PTV and synonymous
variants in 31,058 trios, and fits BurdenMLE-DN with the bundled calibrated
mutation rates and six LOEUF/mutation-rate strata. It writes one results table
and a simple PDF with 95% gene-bootstrap intervals (100 replicates).

Install the package and plotting dependency, then choose an empty working directory:

```r
install.packages(c("remotes", "ggplot2"))
remotes::install_github("ajaynadig/BurdenMLE-DN")
```

If testing an unmerged example branch, install that branch rather than `main`.
From a clone, `R CMD INSTALL .` installs the checked-out package and example.

```r
example_script <- system.file("examples", "kaplanis",
                             "run_kaplanis_public_example.R", package = "BurdenMLEDN")
file.copy(example_script, "run_kaplanis_public_example.R")
```

Open the copied script and run it line by line in RStudio, or use
`source("run_kaplanis_public_example.R")`. Downloads and outputs go into the
working directory, not the installed package. The download overwrites an
existing copy of the public table. Internet access is required; no manuscript
count data are distributed.

The reference contains 17,395 autosomal genes; only variants matching these
gene symbols contribute to the fits. Unobserved reference genes have zero
counts. PTVs comprise the table's frameshift, stop-gained, splice-donor and
splice-acceptor annotations. Reference rates are synonymous-calibrated SNV
rates; an observed PTV/SNV ratio scales PTV rates to include indels.

This public example gives approximately 7.2%
PTV mutational variance versus approximately 6.5% in the main analysis, likely due to annotation/variant-mask and SNV/indel-handling differences.

For your own counts, replace the download/counting section with a data frame
in reference gene order containing `case_count`, per-haploid `case_rate`, and
`N`. Use compatible variant masks and callable exposure; use
`loeuf_features(reference$gene_id, reference)` for the feature matrix.
