# BurdenMLE-DN

BurdenMLE-DN fits a discrete-mixture maximum-likelihood model to gene-level
de novo variant counts. It estimates the distribution of gene-level rate
ratios and derived quantities including observed-scale mutational variance and
effective penetrance.

## Installation and quick start

```r
install.packages("remotes")
remotes::install_github("ajaynadig/BurdenMLE-DN")

library(BurdenMLEDN)
example_file <- system.file(
  "extdata", "synthetic_example.csv", package = "BurdenMLEDN"
)
dat <- read.csv(example_file, row.names = "gene")

fit <- BurdenMLE_DN(
  dat,
  prevalence = 0.02,
  bootstrap = FALSE,
  null_sim = FALSE
)
fit
mixture_weights(fit)
```

The bundled example is entirely synthetic. For uncertainty estimates, set
`bootstrap = TRUE` and choose `n_boot`; null simulations are controlled
separately by `null_sim` and `n_null`.

## Documentation

- [Input data and model parameters](wiki/Input-Data-and-Model-Parameters.md)
- [Output fields](wiki/Output-Fields.md)
- [Worked example](wiki/Worked-Example.md)
- [Reproducing the manuscript](wiki/Reproducing-the-Manuscript.md)
- [Citation](wiki/Citation.md)

The package does **not** distribute individual-level variants, gene-level
counts from the manuscript cohorts, or other controlled/embargoed inputs.
Those inputs remain under the governance of the originating study authors.

## License

MIT © 2026 Ajay Nadig.
