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
features <- model.matrix(~ analysis_stratum - 1, dat)
rownames(features) <- rownames(dat)

fit <- BurdenMLE_DN(
  dat,
  features = features,
  prevalence = 0.02,
  bootstrap = FALSE
)
fit
summary(fit)
```

The bundled data are simulated and contain no study-participant information.
LOEUF-like annotation strata are included because substantive analyses should
generally allow the mixture weights to differ across gene-constraint strata.

## Pages

- [Input data and model parameters](Input-Data-and-Model-Parameters)
- [Output fields](Output-Fields)
- [Worked example](Worked-Example)
- [Reproducing the manuscript](Reproducing-the-Manuscript)
- [Citation](Citation)
