# Worked example

The installed package includes a reproducible synthetic dataset.

```r
library(BurdenMLEDN)

path <- system.file("extdata", "synthetic_example.csv",
                    package = "BurdenMLEDN")
dat <- read.csv(path, row.names = "gene")

set.seed(1)
fit <- BurdenMLE_DN(
  input_data = dat,
  prevalence = 0.02,
  no_cpts = 10,
  grid_size = 10,
  bootstrap = TRUE,
  n_boot = 20,
  null_sim = FALSE
)

fit
fit$mutvar_output$total_mutvar
fit$mutvar_output$mutvar_CI
fit$penetrance
mixture_weights(fit)
```

Twenty bootstrap replicates keep the example short; use at least the
manuscript-specified count for substantive inference. The bootstrap resamples
genes and refits the mixture in each replicate.

For annotation-specific weights:

```r
features <- model.matrix(~ annotation - 1, dat)
rownames(features) <- rownames(dat)

fit_by_annotation <- BurdenMLE_DN(
  dat,
  features = features,
  prevalence = 0.02,
  bootstrap = FALSE,
  null_sim = FALSE
)
mixture_weights(fit_by_annotation)
```
