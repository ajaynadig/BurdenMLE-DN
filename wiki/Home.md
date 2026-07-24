# BurdenMLE-DN

BurdenMLE-DN estimates a mixture distribution of gene-level de novo variant
rate ratios and derived quantities such as mutational variance and effective
penetrance.

## Installation

```r
install.packages("remotes")
remotes::install_github("ajaynadig/BurdenMLE-DN")
```

## Quick start

```r
library(BurdenMLEDN)
path <- system.file("extdata", "synthetic_example.csv",
                    package = "BurdenMLEDN")
dat <- read.csv(path, row.names = "gene")

fit <- BurdenMLE_DN(
  dat,
  prevalence = 0.02,
  bootstrap = FALSE,
  null_sim = FALSE
)
fit
summary(fit)
```

The bundled data are simulated and contain no study-participant information.

## Pages

- [Input data and model parameters](Input-Data-and-Model-Parameters)
- [Output fields](Output-Fields)
- [Worked example](Worked-Example)
- [Reproducing the manuscript](Reproducing-the-Manuscript)
- [Citation](Citation)
