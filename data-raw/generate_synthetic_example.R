# Generate the package's fully synthetic worked example.
set.seed(20260724)

n_genes <- 240L
sample_size <- 15000L
gene <- sprintf("SYNTHETIC_GENE_%03d", seq_len(n_genes))
annotation <- rep(c("PTV_like", "missense_like"), each = n_genes / 2)
case_rate <- exp(rnorm(n_genes, log(1.4e-5), 0.45))
loeuf_quintile <- rep(seq_len(5), length.out = n_genes)
high_mutation_rate <- case_rate > stats::quantile(case_rate, 0.8)
analysis_stratum <- ifelse(
  loeuf_quintile == 1,
  ifelse(high_mutation_rate, "LOEUF1_mu1", "LOEUF1_mu2"),
  paste0("LOEUF", loeuf_quintile)
)

# Simulated rate ratios are used only to create an illustrative dataset.
log_rr <- ifelse(
  runif(n_genes) < 0.12,
  runif(n_genes, log(1.5), log(6)),
  0
)
case_count <- rpois(n_genes, 2 * sample_size * case_rate * exp(log_rr))

synthetic_example <- data.frame(
  gene = gene,
  annotation = annotation,
  loeuf_quintile = loeuf_quintile,
  analysis_stratum = analysis_stratum,
  case_count = case_count,
  case_rate = case_rate,
  N = sample_size,
  check.names = FALSE
)

write.csv(
  synthetic_example,
  file.path("inst", "extdata", "synthetic_example.csv"),
  row.names = FALSE,
  quote = FALSE
)
