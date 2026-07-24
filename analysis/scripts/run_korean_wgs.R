#!/usr/bin/env Rscript

# Reproduce the Kim et al. Korean autism WGS analysis from annotated de novo
# variants and the same mutation-rate references used by the main study.
# Model fitting and bootstrapping use the production BurdenMLE-DN code.

options(stringsAsFactors = FALSE, warn = 1)
library(pbapply)

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_file <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
} else {
  normalizePath(sys.frame(1)$ofile, mustWork = TRUE)
}
example_dir <- dirname(script_file)
repo_dir <- dirname(example_dir)
final_runs_dir <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT", unset = file.path(repo_dir, "analysis")), mustWork = FALSE)

parse_args <- function(x) {
  out <- list(
    n_boot = 100L,
    seed = 20260721L,
    run_date = format(Sys.Date(), "%b%d_%y")
  )
  i <- 1L
  while (i <= length(x)) {
    key <- sub("^--", "", x[[i]])
    if (!key %in% c("n-boot", "seed", "run-date")) {
      stop("Unknown argument: ", x[[i]])
    }
    if (i == length(x)) stop("Missing value after ", x[[i]])
    value <- x[[i + 1L]]
    key <- gsub("-", "_", key)
    out[[key]] <- value
    i <- i + 2L
  }
  out$n_boot <- as.integer(out$n_boot)
  out$seed <- as.integer(out$seed)
  if (is.na(out$n_boot) || out$n_boot < 1L) {
    stop("--n-boot must be a positive integer.")
  }
  if (is.na(out$seed)) stop("--seed must be an integer.")
  if (!grepl("^[A-Za-z0-9_-]+$", out$run_date)) {
    stop("--run-date may contain only letters, numbers, underscores, and hyphens.")
  }
  out
}

opt <- parse_args(args)

input_dir <- file.path(final_runs_dir, "inputs")
derived_dir <- file.path(final_runs_dir, "outputs", "derived", "korean_wgs")
model_dir <- file.path(
  final_runs_dir, "outputs", "models", "sensitivities", "korean_wgs"
)
figure_dir <- file.path(final_runs_dir, "outputs", "figures", "supplementary")
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Make every stochastic step reproducible, including the gene bootstrap.
set.seed(opt$seed)

required_files <- c(
  kim_variants = file.path(input_dir, "korean_wgs", "kim_2024.tsv"),
  mutation_rates = file.path(
    input_dir, "mutation_rates",
    "ASD_gene_table_w_bespoke_mutation_rates_2024-07-24.txt"
  ),
  gnomad_constraint = file.path(
    input_dir, "constraint", "gnomad.v4.1.constraint_metrics.tsv"
  )
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required source file(s):\n", paste(" -", missing_files, collapse = "\n"))
}
if (!requireNamespace("mixsqp", quietly = TRUE)) {
  stop("The R package 'mixsqp' is required. Install it with install.packages('mixsqp').")
}

source_files <- c(
  "BurdenMLE_DN.R", "EM.R", "MixSQP.R", "estimate_mutvar.R", "io.R",
  "likelihoods.R", "model.R", "secondary_analysis_functions.R"
)
for (source_file in source_files) {
  source(if (identical(source_file, "secondary_analysis_functions.R")) file.path(repo_dir, "analysis", "scripts", source_file) else file.path(repo_dir, "R", source_file))
}

message("Kim et al. (2024) Korean autism WGS analysis")
message("Started: ", format(Sys.time(), tz = "UTC", usetz = TRUE))

read_tsv <- function(path, select = NULL) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    out <- data.table::fread(path, select = select, data.table = FALSE)
  } else {
    out <- utils::read.delim(path, check.names = FALSE, quote = "", comment.char = "")
    if (!is.null(select)) out <- out[, select, drop = FALSE]
  }
  out
}

resolve_duplicate_gene_ids <- function(x) {
  replacements <- data.frame(
    Gene = c("HIST1H4F", "ATXN7", "PRSS50", "CCDC39", "HSPA14", "IGF2", "PDE11A"),
    Gene_ID = c("ENSG00000274618", "ENSG00000163635", "ENSG00000283706",
                "ENSG00000284862", "ENSG00000187522", "ENSG00000167244",
                "ENSG00000128655"),
    LOEUF = c(NA, 0.328, 0.99, 0.754, 0.395, 1.131, 1.404)
  )
  for (i in seq_len(nrow(replacements))) {
    idx <- x$Gene == replacements$Gene[[i]]
    x$Gene_ID[idx] <- replacements$Gene_ID[[i]]
    x$LOEUF[idx] <- replacements$LOEUF[[i]]
  }
  multi_chr <- x$Gene %in% c("RF00017", "RF00019")
  x$Gene[multi_chr] <- paste0(x$Gene[multi_chr], "_chr", x$Chrom[multi_chr])
  rf19 <- c(
    RF00019_chr3 = "ENSG00000212392", RF00019_chr6 = "ENSG00000200314",
    RF00019_chr7 = "ENSG00000201913", RF00019_chr12 = "ENSG00000207176",
    RF00019_chr16 = "ENSG00000199668"
  )
  for (nm in names(rf19)) {
    idx <- x$Gene == nm
    x$Gene_ID[idx] <- rf19[[nm]]
    x$LOEUF[idx] <- NA_real_
  }
  x
}

count_by_gene <- function(gene_ids, keep, target_ids) {
  tabulate(match(gene_ids[keep], target_ids), nbins = length(target_ids))
}

message("Reading annotated de novo variants...")
kim <- read_tsv(required_files[["kim_variants"]])
required_kim_columns <- c(
  "chr", "variant_class", "vcf_iid", "vcf_fid",
  "most_severe_consequence", "gene_symbol", "gene_id", "loeuf"
)
missing_kim_columns <- setdiff(required_kim_columns, names(kim))
if (length(missing_kim_columns)) {
  stop(
    "The Korean WGS variant table is missing required columns: ",
    paste(missing_kim_columns, collapse = ", ")
  )
}

# The published supplementary variant table includes probands and other
# sequenced children. In these family IDs, the child-role suffix immediately
# before "-blood" is 3 for the designated proband and 4-6 for additional
# children. The 680 suffix-3 individuals are the proband set described in the
# reviewer response.
source_variant_rows <- nrow(kim)
sample_core_id <- sub("-blood.*$", "", kim$vcf_iid)
child_role_suffix <- substring(sample_core_id, nchar(sample_core_id))
sample_role_counts <- table(
  child_role_suffix[!duplicated(kim$vcf_iid)]
)
kim <- kim[child_role_suffix == "3", , drop = FALSE]

ptv_consequences <- c(
  "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant",
  "stop_gained", "stop_lost", "start_lost"
)
kim$isPTV <- kim$most_severe_consequence %in% ptv_consequences
kim$isMis <- kim$most_severe_consequence == "missense_variant"
kim$isSyn <- kim$most_severe_consequence == "synonymous_variant"
kim$isIndel <- kim$variant_class == "INDEL"
kim$Gene <- kim$gene_symbol
kim$Gene_ID <- kim$gene_id
kim$LOEUF <- suppressWarnings(as.numeric(kim$loeuf))
kim$Chrom <- sub("^chr", "", kim$chr)
kim <- resolve_duplicate_gene_ids(kim)

n_probands <- length(unique(kim$vcf_iid))
if (n_probands != length(unique(kim$vcf_fid))) {
  stop("Proband individual and family counts do not agree after role filtering.")
}
ptv_total <- sum(kim$isPTV, na.rm = TRUE)
ptv_snv_total <- sum(kim$isPTV & !kim$isIndel, na.rm = TRUE)
ptv_scale <- ptv_total / ptv_snv_total
message(
  "Source variants: ", source_variant_rows,
  "; proband variants after role filtering: ", nrow(kim),
  "; unique affected offspring: ", n_probands
)
message("PTVs: ", ptv_total, " total; ", ptv_snv_total, " SNVs; scale factor: ", signif(ptv_scale, 5))

message("Reading mutation rates and constructing gene table...")
mut <- read_tsv(required_files[["mutation_rates"]])
rate_cols <- setdiff(names(mut), names(mut)[1:5])
gene <- data.frame(
  gene_id = mut$Gene_ID,
  gene = mut$Gene,
  LOEUF = suppressWarnings(as.numeric(mut$LOEUF)),
  Chrom = as.character(mut$Chrom),
  mut[, rate_cols, drop = FALSE],
  check.names = FALSE
)

# The gene universe must come from the mutation-rate table so that genes with
# zero Korean de novo variants remain in the likelihood. Map variants by stable
# Ensembl ID, with gene symbol as a fallback for legacy ID mismatches.
kim_gene_index <- match(kim$Gene_ID, gene$gene_id)
fallback <- is.na(kim_gene_index) & !is.na(kim$Gene) & nzchar(kim$Gene)
kim_gene_index[fallback] <- match(kim$Gene[fallback], gene$gene)
kim$AnalysisGeneID <- gene$gene_id[kim_gene_index]
valid_mapping <- !is.na(kim$AnalysisGeneID)

gene$PTV <- count_by_gene(kim$AnalysisGeneID, kim$isPTV & !kim$isIndel & valid_mapping, gene$gene_id)
gene$PTV_All <- count_by_gene(kim$AnalysisGeneID, kim$isPTV & valid_mapping, gene$gene_id)
gene$MisAll <- count_by_gene(
  kim$AnalysisGeneID, kim$isMis & valid_mapping, gene$gene_id
)
gene$Syn <- count_by_gene(kim$AnalysisGeneID, kim$isSyn & valid_mapping, gene$gene_id)

gene <- gene[complete.cases(gene[, rate_cols, drop = FALSE]), , drop = FALSE]
gene <- gene[!is.na(gene$LOEUF), , drop = FALSE]
gene <- gene[!gene$Chrom %in% c("X", "Y"), , drop = FALSE]
rownames(gene) <- gene$gene_id
message("Genes before gnomAD calibration match: ", nrow(gene))

q <- stats::quantile(gene$LOEUF, probs = seq(0, 1, 0.2))
loeuf_bin <- as.integer(cut(gene$LOEUF, breaks = q, include.lowest = TRUE, labels = 1:5))
names(loeuf_bin) <- gene$gene_id

message("Calibrating gene-wise mutation rates with gnomAD v4.1 synonymous counts...")
gnomad <- read_tsv(
  required_files[["gnomad_constraint"]],
  select = c("gene_id", "canonical", "syn.obs", "syn.exp")
)
gkeep <- gnomad$canonical %in% TRUE & grepl("ENSG", gnomad$gene_id) &
  complete.cases(gnomad[, c("gene_id", "syn.obs", "syn.exp")])
gdat <- data.frame(
  case_count = gnomad$syn.obs[gkeep],
  expected_count = gnomad$syn.exp[gkeep],
  gene_id = gnomad$gene_id[gkeep]
)
rownames(gdat) <- gnomad$gene_id[gkeep]
cal_model <- BurdenMLE_DN(
  input_data = gdat,
  features = NULL,
  component_endpoints = seq(-2, 2, length.out = 31),
  mutvar_est = FALSE,
  null_sim = FALSE,
  bootstrap = FALSE,
  return_likelihood = TRUE,
  estimate_posteriors = TRUE,
  estimate_effective_penetrance = FALSE,
  optimizer = "mixsqp"
)
correction <- cal_model$posterior_gene_estimates$Posterior_Mean
names(correction) <- rownames(cal_model$posterior_gene_estimates)
gene$PosteriorMuCorrectionFactor <- correction[gene$gene_id]
gene <- gene[!is.na(gene$PosteriorMuCorrectionFactor), , drop = FALSE]
loeuf_bin <- loeuf_bin[match(gene$gene_id, names(loeuf_bin))]
message("Final analyzed genes: ", nrow(gene))

mu_original <- gene[, c("mu_snp_PTV", "mu_snp_Mis2", "mu_snp_Mis1", "mu_snp_Mis0", "mu_snp_Syn")]
cumulative_mu <- rowSums(mu_original) * gene$PosteriorMuCorrectionFactor
high_mu <- cumulative_mu > stats::quantile(cumulative_mu, 0.8)
features <- cbind(
  LOEUF1_mu1 = as.numeric(loeuf_bin == 1 & high_mu),
  LOEUF1_mu2 = as.numeric(loeuf_bin == 1 & !high_mu),
  LOEUF2 = as.numeric(loeuf_bin == 2),
  LOEUF3 = as.numeric(loeuf_bin == 3),
  LOEUF4 = as.numeric(loeuf_bin == 4),
  LOEUF5 = as.numeric(loeuf_bin == 5)
)
rownames(features) <- gene$gene_id
stopifnot(all(rowSums(features) == 1))
gene$Stratum <- sub("^LOEUF", "OE", colnames(features)[max.col(features)])

# Pool the three missense mutation-rate categories because AlphaMissense
# annotations are unavailable in the published Korean variant table.
gene$mu_snp_MisAll <- gene$mu_snp_Mis2 + gene$mu_snp_Mis1 + gene$mu_snp_Mis0
syn_expected <- sum(
  2 * n_probands * gene$mu_snp_Syn * gene$PosteriorMuCorrectionFactor
)
message(
  "Synonymous negative-control O/E: ",
  signif(sum(gene$Syn) / syn_expected, 7)
)

component_endpoints <- seq(0, log(100), length.out = 10)
prevalence <- 0.0276
classes <- data.frame(
  class = c("PTV", "Pooled missense", "Synonymous"),
  count_col = c("PTV", "MisAll", "Syn"),
  rate_col = c("mu_snp_PTV", "mu_snp_MisAll", "mu_snp_Syn"),
  stringsAsFactors = FALSE
)

make_class_input <- function(count_col, rate_col) {
  case_rate <- gene[[rate_col]] * gene$PosteriorMuCorrectionFactor
  keep <- !(gene[[count_col]] > 0 & case_rate == 0)
  if (!all(keep)) {
    stop(
      "Observed ", count_col,
      " variants were assigned to a gene with zero mutation rate."
    )
  }
  dat <- data.frame(
    case_count = gene[[count_col]],
    case_rate = case_rate,
    N = n_probands
  )
  rownames(dat) <- gene$gene_id
  dat
}

fit_inputs <- lapply(seq_len(nrow(classes)), function(i) {
  make_class_input(classes$count_col[[i]], classes$rate_col[[i]])
})
names(fit_inputs) <- classes$class

set.seed(opt$seed)
bootstrap_samples <- replicate(
  opt$n_boot, sample.int(nrow(gene), replace = TRUE)
)
message(
  "Fitting PTV, pooled-missense, and synonymous BurdenMLE-DN models with ",
  opt$n_boot, " shared gene-bootstrap replicates..."
)
fits <- lapply(seq_len(nrow(classes)), function(i) {
  message("  ", classes$class[[i]])
  BurdenMLE_DN(
    input_data = fit_inputs[[i]],
    features = features,
    component_endpoints = component_endpoints,
    prevalence = prevalence,
    mutvar_est = TRUE,
    bootstrap = TRUE,
    bootstrap_samples = bootstrap_samples,
    n_boot = opt$n_boot,
    null_sim = FALSE,
    estimate_posteriors = FALSE,
    estimate_effective_penetrance = FALSE,
    optimizer = "mixsqp"
  )
})
names(fits) <- classes$class

all_converged <- vapply(fits, function(model) {
  full_converged <- all(vapply(
    model$mixsqp_output,
    function(x) identical(x$status, "converged to optimal solution"),
    logical(1)
  ))
  full_converged && all(vapply(
    model$bootstrap_output$bootstrap_delta,
    function(x) all(is.finite(x)),
    logical(1)
  ))
}, logical(1))
if (!all(all_converged)) {
  stop(
    "At least one Korean WGS model or bootstrap fit failed: ",
    paste(names(all_converged)[!all_converged], collapse = ", ")
  )
}

likelihood_observed <- vapply(
  fit_inputs, function(x) sum(x$case_count), numeric(1)
)
likelihood_expected <- vapply(
  fit_inputs, function(x) sum(2 * x$N * x$case_rate), numeric(1)
)
observed <- likelihood_observed
expected <- likelihood_expected
observed[[1]] <- sum(gene$PTV_All)
expected[[1]] <- likelihood_expected[[1]] * ptv_scale
oe <- observed / expected
oe_ci <- t(vapply(seq_along(observed), function(i) {
  stats::poisson.test(observed[[i]], T = expected[[i]])$conf.int
}, numeric(2)))

reported_scale <- c(ptv_scale, 1, 1)
point_unscaled <- vapply(
  fits, function(model) model$mutvar_output$total_mutvar, numeric(1)
)
point <- point_unscaled * reported_scale
ci <- t(vapply(seq_along(fits), function(i) {
  fits[[i]]$mutvar_output$mutvar_CI * reported_scale[[i]]
}, numeric(2)))
results <- data.frame(
  variant_class = classes$class,
  n_probands = n_probands,
  n_genes = nrow(gene),
  likelihood_observed = likelihood_observed,
  likelihood_expected = likelihood_expected,
  observed = observed,
  expected = expected,
  observed_expected = oe,
  oe_ci_lower = oe_ci[, 1],
  oe_ci_upper = oe_ci[, 2],
  mutational_variance = point,
  mutvar_ci_lower = ci[, 1],
  mutvar_ci_upper = ci[, 2],
  ptv_indel_scale = reported_scale,
  prevalence = prevalence,
  log_likelihood = vapply(fits, function(x) x$ll, numeric(1)),
  active_components_min = vapply(
    fits, function(x) min(rowSums(x$delta > 1e-8)), numeric(1)
  ),
  active_components_mean = vapply(
    fits, function(x) mean(rowSums(x$delta > 1e-8)), numeric(1)
  ),
  active_components_max = vapply(
    fits, function(x) max(rowSums(x$delta > 1e-8)), numeric(1)
  ),
  n_boot = opt$n_boot,
  seed = opt$seed,
  optimizer = "mixsqp",
  stringsAsFactors = FALSE
)

results_file <- file.path(derived_dir, "korean_wgs_results.tsv")
utils::write.table(
  results, results_file, sep = "\t", quote = FALSE, row.names = FALSE
)

input_summary <- data.frame(
  metric = c(
    "source_variant_rows_all_children",
    "source_unique_proband_individuals",
    "source_unique_nonproband_children",
    "proband_variant_rows",
    "proband_PTV_total",
    "proband_PTV_SNV",
    "proband_pooled_missense",
    "proband_synonymous",
    "proband_coding_variants_mapped_to_reference",
    "analyzed_autosomal_genes",
    "analyzed_PTV_total",
    "analyzed_PTV_SNV",
    "analyzed_pooled_missense",
    "analyzed_synonymous",
    "synonymous_observed_negative_control",
    "synonymous_expected_negative_control",
    "synonymous_OE_negative_control",
    paste0("genes_", colnames(features))
  ),
  value = c(
    source_variant_rows,
    unname(sample_role_counts[["3"]]),
    sum(sample_role_counts[names(sample_role_counts) != "3"]),
    nrow(kim),
    ptv_total,
    ptv_snv_total,
    sum(kim$isMis),
    sum(kim$isSyn),
    sum(
      valid_mapping & (kim$isPTV | kim$isMis | kim$isSyn),
      na.rm = TRUE
    ),
    nrow(gene),
    sum(gene$PTV_All),
    sum(gene$PTV),
    sum(gene$MisAll),
    sum(gene$Syn),
    sum(gene$Syn),
    syn_expected,
    sum(gene$Syn) / syn_expected,
    colSums(features)
  ),
  stringsAsFactors = FALSE
)
input_summary_file <- file.path(
  derived_dir, "korean_wgs_input_summary.tsv"
)
utils::write.table(
  input_summary, input_summary_file,
  sep = "\t", quote = FALSE, row.names = FALSE
)

model_file <- file.path(
  model_dir,
  paste0("korean_wgs_models_mixsqp_", opt$run_date, ".rds")
)
saveRDS(
  list(
    models = fits,
    mutation_rate_calibration_model = cal_model,
    inputs = fit_inputs,
    features = features,
    results = results,
    input_summary = input_summary,
    component_endpoints = component_endpoints,
    prevalence = prevalence,
    bootstrap_samples = bootstrap_samples,
    run_date = opt$run_date,
    seed = opt$seed
  ),
  model_file
)

plot_figure <- function(device, filename) {
  device(filename, width = 7.2, height = 3.35)
  on.exit(grDevices::dev.off())
  layout(matrix(1:2, nrow = 1), widths = c(0.92, 1.08))
  cols <- c("#FF5F8A", "#F9B332", "grey50")
  labels <- c("PTV", "Missense", "Synonymous")

  par(family = "Helvetica", mar = c(4.6, 4.3, 1.1, 0.8), mgp = c(2.35, 0.65, 0),
      tcl = -0.22, las = 1)
  yr <- range(c(oe_ci, 1), finite = TRUE)
  pad <- diff(yr) * 0.12
  plot(seq_along(oe), oe, type = "n", xaxt = "n", xlab = "", ylab = "Observed / expected\n(95% CI)",
       ylim = c(max(0, yr[1] - pad), yr[2] + pad), xlim = c(0.55, 3.45), bty = "l")
  abline(h = 1, col = "#888888", lty = 2, lwd = 1)
  segments(seq_along(oe), oe_ci[, 1], seq_along(oe), oe_ci[, 2], col = cols, lwd = 1.6)
  segments(seq_along(oe) - 0.06, oe_ci[, 1], seq_along(oe) + 0.06, oe_ci[, 1], col = cols, lwd = 1.6)
  segments(seq_along(oe) - 0.06, oe_ci[, 2], seq_along(oe) + 0.06, oe_ci[, 2], col = cols, lwd = 1.6)
  points(seq_along(oe), oe, pch = 21, bg = cols, col = "white", cex = 1.55, lwd = 0.8)
  axis(1, at = seq_along(oe), labels = labels, tick = FALSE, line = 0.15,
       cex.axis = 0.78, gap.axis = -1)
  mtext("A", side = 3, adj = 0, line = 0.05, font = 2, cex = 1.05)

  par(family = "Helvetica", mar = c(4.6, 4.5, 1.1, 0.6), mgp = c(2.55, 0.65, 0),
      tcl = -0.22, las = 1)
  top <- max(ci[, 2]) * 1.13
  plot(seq_along(point), point, type = "n", xaxt = "n", xlab = "", ylab = "Mutational variance\n(95% CI)",
       ylim = c(0, top), xlim = c(0.55, 3.45), bty = "l")
  abline(h = 0, col = "#888888", lty = 2, lwd = 1)
  segments(seq_along(point), ci[, 1], seq_along(point), ci[, 2], col = cols, lwd = 1.6)
  segments(seq_along(point) - 0.06, ci[, 1], seq_along(point) + 0.06, ci[, 1], col = cols, lwd = 1.6)
  segments(seq_along(point) - 0.06, ci[, 2], seq_along(point) + 0.06, ci[, 2], col = cols, lwd = 1.6)
  points(seq_along(point), point, pch = 21, bg = cols, col = "white", cex = 1.55, lwd = 0.8)
  axis(1, at = seq_along(point), labels = labels, tick = FALSE, line = 0.15,
       cex.axis = 0.78, gap.axis = -1)
  mtext("B", side = 3, adj = 0, line = 0.05, font = 2, cex = 1.05)
}

figure_pdf <- file.path(
  figure_dir, "SupplementaryFigure7_KoreanWGS.pdf"
)
plot_figure(grDevices::pdf, figure_pdf)
message("\nFinal results:")
print(results, digits = 5, row.names = FALSE)
message("Results: ", results_file)
message("Input summary: ", input_summary_file)
message("Models: ", model_file)
message("Figure: ", figure_pdf)
message("Completed: ", format(Sys.time(), tz = "UTC", usetz = TRUE))
