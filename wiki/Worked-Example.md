# Worked example

The installed package includes a reproducible synthetic dataset.

```r
library(BurdenMLEDN)

path <- system.file("extdata", "synthetic_example.csv",
                    package = "BurdenMLEDN")
dat <- read.csv(path, row.names = "gene")
features <- model.matrix(~ analysis_stratum - 1, dat)
rownames(features) <- rownames(dat)

set.seed(1)
fit <- BurdenMLE_DN(
  input_data = dat,
  features = features,
  prevalence = 0.02,
  no_cpts = 10,
  grid_size = 10,
  bootstrap = TRUE,
  n_boot = 20
)

fit
fit$mutvar_output$total_mutvar
fit$mutvar_output$mutvar_CI
fit$penetrance
mixture_weights(fit)

sampler <- posterior_gene_sampler(fit, dat)
log_rate_ratios <- posterior_gene_samples(sampler, num_samples = 3)
```

Twenty bootstrap replicates keep the example short; use at least the
manuscript-specified count for substantive inference. The bootstrap resamples
genes and refits the mixture in each replicate. Set `bootstrap_seed` to replay
generated samples, or supply a character `bootstrap_samples` matrix containing
gene row names when several related fits must use exactly the same draws.

The prepared posterior sampler is intended for repeated simulations such as
gene-discovery forecasts. It retains the fitted component probabilities, then
conditions each continuous within-component draw on that gene's observed and
expected counts. Its result is always a named gene-by-sample matrix. Set the
ordinary R seed before sampling when draws must be reproducible.

For a real analysis, the bundled count-free reference provides the baseline
LOEUF annotations and mutation rates:

```r
reference <- burdenmle_gene_reference()
features <- loeuf_features(reference$gene_id)

# After matching your own observed counts into this gene order:
input <- data.frame(
  case_count = your_observed_counts,
  case_rate = reference$corrected_mu_snp_Mis2,
  N = rep(your_number_of_trios, nrow(reference)),
  row.names = reference$gene_id
)

fit <- BurdenMLE_DN(input, features = features, prevalence = 0.02)
```

`your_observed_counts` and `your_number_of_trios` are placeholders for data
supplied by the user; no manuscript counts are bundled with the package.
