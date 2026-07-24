# Export only non-sensitive gene annotations and mutation-rate resources from
# a completed manuscript model file. No assays, cohort counts, sample metadata,
# or fitted model objects are written.

args <- commandArgs(trailingOnly = TRUE)
model_file <- if (length(args)) args[1] else Sys.getenv(
  "BURDENMLEDN_MODEL_FILE", unset = ""
)
if (!nzchar(model_file) || !file.exists(model_file)) {
  stop(
    "Supply the authorized main autism model file as the first argument or ",
    "BURDENMLEDN_MODEL_FILE."
  )
}

environment <- new.env(parent = emptyenv())
load(model_file, envir = environment)
if (!exists("autism_data", envir = environment, inherits = FALSE)) {
  stop("The model file does not contain autism_data.")
}

data <- environment$autism_data
if (!methods::is(data$counts, "SummarizedExperiment")) {
  stop("autism_data$counts is not a SummarizedExperiment.")
}
row_data <- as.data.frame(SummarizedExperiment::rowData(data$counts))
features <- as.matrix(data$features)
if (!identical(rownames(row_data), rownames(features))) {
  stop("Gene annotations and features are not aligned.")
}
if (!all(rowSums(features) == 1)) {
  stop("The manuscript features are not one-hot.")
}

loeuf_columns <- paste0("LOEUF", 1:5)
rate_columns <- c(
  "mu_snp_PTV", "mu_snp_Mis2", "mu_snp_Mis1",
  "mu_snp_Mis0", "mu_snp_Syn"
)
required <- c(
  "gene_id", "gene", "Chrom", "LOEUF",
  loeuf_columns, rate_columns, "PosteriorMuCorrectionFactor"
)
missing_columns <- setdiff(required, names(row_data))
if (length(missing_columns)) {
  stop(
    "Required public annotation column(s) are absent: ",
    paste(missing_columns, collapse = ", ")
  )
}

loeuf_quintile <- max.col(as.matrix(row_data[, loeuf_columns]))
analysis_stratum <- colnames(features)[max.col(features)]
correction <- row_data$PosteriorMuCorrectionFactor

reference <- data.frame(
  gene_id = row_data$gene_id,
  gene_symbol = row_data$gene,
  chromosome = row_data$Chrom,
  loeuf = row_data$LOEUF,
  loeuf_quintile = loeuf_quintile,
  analysis_stratum = analysis_stratum,
  mutation_rate_correction = correction,
  row_data[, rate_columns, drop = FALSE],
  check.names = FALSE
)
for (column in rate_columns) {
  reference[[paste0("corrected_", column)]] <-
    reference[[column]] * reference$mutation_rate_correction
}
reference <- reference[order(reference$gene_id), ]

allowed_strata <- c(
  "LOEUF1_mu1", "LOEUF1_mu2",
  "LOEUF2", "LOEUF3", "LOEUF4", "LOEUF5"
)
if (!setequal(unique(reference$analysis_stratum), allowed_strata)) {
  stop("Unexpected analysis stratum encountered.")
}
if (anyNA(reference) || anyDuplicated(reference$gene_id)) {
  stop("The public gene reference is incomplete or has duplicate gene IDs.")
}

dir.create(file.path("inst", "extdata"), recursive = TRUE, showWarnings = FALSE)
write.csv(
  reference,
  file.path("inst", "extdata", "burdenmle_gene_reference.csv"),
  row.names = FALSE,
  quote = FALSE
)

description <- c(
  LOEUF1_mu1 = paste(
    "Most constrained LOEUF quintile and above the global 80th percentile",
    "of corrected cumulative coding SNV mutation rate"
  ),
  LOEUF1_mu2 = paste(
    "Most constrained LOEUF quintile and at or below the global 80th",
    "percentile of corrected cumulative coding SNV mutation rate"
  ),
  LOEUF2 = "Second LOEUF quintile",
  LOEUF3 = "Third LOEUF quintile",
  LOEUF4 = "Fourth LOEUF quintile",
  LOEUF5 = "Least constrained LOEUF quintile"
)
strata <- do.call(rbind, lapply(allowed_strata, function(stratum) {
  selected <- reference$analysis_stratum == stratum
  data.frame(
    analysis_stratum = stratum,
    description = unname(description[stratum]),
    genes = sum(selected),
    loeuf_min = min(reference$loeuf[selected]),
    loeuf_max = max(reference$loeuf[selected]),
    stringsAsFactors = FALSE
  )
}))
write.csv(
  strata,
  file.path("inst", "extdata", "loeuf_strata_summary.csv"),
  row.names = FALSE,
  quote = TRUE
)
