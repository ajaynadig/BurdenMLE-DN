# Build an ordered, manuscript-facing checklist of study-derived numerical
# estimates. This script reads finalized derived outputs and fitted models; it
# does not edit the manuscript or refit any model.

suppressPackageStartupMessages(library(SummarizedExperiment))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit)) return(sub(paste0("^", flag, "="), "", hit[1]))
  index <- match(flag, args)
  if (!is.na(index)) return(args[index + 1])
  default
}

script_file <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
)
script_file <- normalizePath(script_file, mustWork = TRUE)
repo_dir <- dirname(dirname(dirname(script_file)))
final_runs_dir <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT", unset = file.path(repo_dir, "analysis")), mustWork = FALSE)
outputs_dir <- file.path(final_runs_dir, "outputs")
derived_dir <- file.path(outputs_dir, "derived")
tables_dir <- file.path(outputs_dir, "tables")
output_dir <- get_arg(
  "--output-dir",
  file.path(tables_dir, "manuscript_values")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Required derived output is missing: ", path)
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}
read_csv <- function(path) {
  if (!file.exists(path)) stop("Required derived output is missing: ", path)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}
relative_path <- function(path) {
  normalized <- normalizePath(path, mustWork = FALSE)
  prefix <- paste0(normalizePath(final_runs_dir), "/")
  sub(paste0("^", prefix), "", normalized)
}
one_row <- function(data, ...) {
  result <- subset(data, ...)
  if (nrow(result) != 1L) {
    stop("Expected one matching row; found ", nrow(result), ".")
  }
  result
}

main_dir <- file.path(derived_dir, "main_autism")
autism_dd_dir <- file.path(derived_dir, "autism_dd")
cohort_ddid_dir <- file.path(derived_dir, "cohort_ddid")
korean_dir <- file.path(derived_dir, "korean_wgs")
paternal_dir <- file.path(tables_dir, "supplementary", "paternal_age")
simulation_diagnostic_file <- file.path(
  tables_dir, "diagnostics", "simulations", "SimulationOptimizerComparison.tsv"
)

main_mutvar_file <- file.path(main_dir, "mutational_variance.tsv")
main_penetrance_file <- file.path(main_dir, "effective_penetrance.tsv")
main_fraction_file <- file.path(main_dir, "fraction_cases_by_rr.tsv")
main_prevalence_file <- file.path(main_dir, "prevalence_mutational_variance.tsv")
main_cohort_file <- file.path(main_dir, "figure2A_cohort_table.tsv")
tada_file <- file.path(main_dir, "tada_comparison_statistics.tsv")
autism_dd_mutvar_file <- file.path(autism_dd_dir, "mutational_variance.tsv")
autism_dd_rr_file <- file.path(autism_dd_dir, "effective_rr.tsv")
autism_dd_fraction_file <- file.path(autism_dd_dir, "fraction_cases.tsv")
autism_dd_genes_file <- file.path(autism_dd_dir, "number_genes.tsv")
korean_results_file <- file.path(korean_dir, "korean_wgs_results.tsv")
korean_input_file <- file.path(korean_dir, "korean_wgs_input_summary.tsv")
gene_table_file <- file.path(
  tables_dir, "supplementary",
  "SupplementaryTable4_GeneLevelMutationalVariance.tsv"
)
cohort_subgroup_file <- file.path(
  tables_dir, "supplementary",
  "SupplementaryTable2_ConstrainedPTVBurdenByStratum.tsv"
)
paternal_summary_file <- file.path(
  paternal_dir, "spark_vs_us_father_age_summary.csv"
)
paternal_test_file <- file.path(
  paternal_dir, "spark_vs_us_father_age_mean_test.csv"
)
paternal_bias_file <- file.path(
  paternal_dir, "spark_paternal_age_mutvar_sensitivity.csv"
)

main_mutvar <- read_tsv(main_mutvar_file)
main_penetrance <- read_tsv(main_penetrance_file)
main_fraction <- read_tsv(main_fraction_file)
main_prevalence <- read_tsv(main_prevalence_file)
main_cohort <- read_tsv(main_cohort_file)
tada <- read_tsv(tada_file)
autism_dd_mutvar <- read_tsv(autism_dd_mutvar_file)
autism_dd_rr <- read_tsv(autism_dd_rr_file)
autism_dd_fraction <- read_tsv(autism_dd_fraction_file)
autism_dd_genes <- read_tsv(autism_dd_genes_file)
korean_results <- read_tsv(korean_results_file)
korean_input <- read_tsv(korean_input_file)
gene_table <- read_tsv(gene_table_file)
cohort_subgroups <- read_tsv(cohort_subgroup_file)
paternal_summary <- read_csv(paternal_summary_file)
paternal_test <- read_csv(paternal_test_file)
paternal_bias <- read_csv(paternal_bias_file)
optimizer_comparison <- read_tsv(simulation_diagnostic_file)

main_summary_file <- file.path(main_dir, "figure2_summary.rds")
if (!file.exists(main_summary_file)) stop("Required summary is missing: ", main_summary_file)
main_summary <- readRDS(main_summary_file)
posterior_summary_file <- file.path(
  main_dir, "posterior_diagnostics_summary.rds"
)
if (!file.exists(posterior_summary_file)) {
  stop("Required summary is missing: ", posterior_summary_file)
}
posterior_summary <- readRDS(posterior_summary_file)

source(file.path(repo_dir, "R", "io.R"))
source(file.path(repo_dir, "analysis", "scripts", "secondary_analysis_functions.R"))
source(file.path(repo_dir, "analysis", "scripts", "model_artifacts.R"))
model_manifest <- get_arg("--model-manifest", NA_character_)
legacy_autism_file <- get_arg("--legacy-autism-model-file", NA_character_)
legacy_ddd_file <- get_arg("--legacy-ddd-model-file", NA_character_)
autism_env <- new.env(parent = globalenv())
ddd_env <- new.env(parent = globalenv())
autism_info <- load_model_artifact(
  model_manifest, "autism", envir = autism_env,
  legacy_path = legacy_autism_file
)
ddd_info <- load_model_artifact(
  model_manifest, "ddd", envir = ddd_env,
  legacy_path = legacy_ddd_file
)
autism_model_file <- autism_info$path
ddd_model_file <- ddd_info$path

# E[(exp(U)-1)^2] for U ~ Uniform(0, endpoint). This is the deterministic
# version of the mutational-variance moment already used in the prevalence
# sensitivity analysis.
effect_moment <- function(endpoints) {
  result <- numeric(length(endpoints))
  nonzero <- abs(endpoints) > 1e-10
  endpoint <- endpoints[nonzero]
  result[nonzero] <-
    expm1(2 * endpoint) / (2 * endpoint) -
    2 * expm1(endpoint) / endpoint + 1
  result
}
annotation_mutvar_constants <- function(model, genetic_data) {
  moments <- effect_moment(model$component_endpoints)
  vapply(seq_len(ncol(model$features)), function(annotation) {
    rows <- model$features[, annotation] == 1
    sum(rows) * 2 * mean(genetic_data$case_rate[rows]) *
      sum(pmax(model$delta[annotation, ], 0) * moments)
  }, numeric(1))
}
mutvar_constant <- function(model, genetic_data) {
  sum(annotation_mutvar_constants(model, genetic_data))
}
bootstrap_constants <- function(model, genetic_data) {
  vapply(seq_along(model$bootstrap_output$bootstrap_delta), function(iteration) {
    indices <- model$bootstrap_output$bootstrap_indices[, iteration]
    boot_model <- model
    boot_model$features <- model$features[indices, , drop = FALSE]
    boot_model$delta <- model$bootstrap_output$bootstrap_delta[[iteration]]
    mutvar_constant(
      boot_model,
      genetic_data[indices, , drop = FALSE]
    )
  }, numeric(1))
}
paired_combined_mutvar <- function(
    data, fitted_models, prevalence_index, ptv_name, mis2_name, label) {
  model_index <- setNames(seq_along(data$loop_vars$names), data$loop_vars$names)
  ptv_index <- model_index[[ptv_name]]
  mis2_index <- model_index[[mis2_name]]
  models <- fitted_models[[prevalence_index]]
  ptv_model <- models[[ptv_index]]
  mis2_model <- models[[mis2_index]]
  ptv_data <- get_genetic_data(ptv_index, data)$genetic_data
  mis2_data <- get_genetic_data(mis2_index, data)$genetic_data
  if (!identical(
    ptv_model$bootstrap_output$bootstrap_samples,
    mis2_model$bootstrap_output$bootstrap_samples
  )) {
    stop(label, ": PTV and Mis2 bootstrap samples are not paired.")
  }
  prevalence <- data$baseprev * data$prev_factors[prevalence_index]
  prevalence_scale <- prevalence / (1 - prevalence)
  ptv_boot <- bootstrap_constants(ptv_model, ptv_data)
  mis2_boot <- bootstrap_constants(mis2_model, mis2_data)
  n_boot <- min(length(ptv_boot), length(mis2_boot))
  combined_boot <- (
    ptv_boot[seq_len(n_boot)] * data$ptv_scale_factor +
      mis2_boot[seq_len(n_boot)]
  ) * prevalence_scale
  data.frame(
    label = label,
    prevalence = prevalence,
    estimate =
      ptv_model$mutvar_output$total_mutvar * data$ptv_scale_factor +
      mis2_model$mutvar_output$total_mutvar,
    lower = unname(quantile(combined_boot, 0.025)),
    upper = unname(quantile(combined_boot, 0.975)),
    n_boot = n_boot
  )
}

autism_primary_combined <- paired_combined_mutvar(
  autism_env$autism_data,
  autism_env$BurdenMLE_DN_models_autism,
  1L,
  "Combined Probands PTV",
  "Combined Probands Mis2",
  "Autism, primary prevalence"
)
one_percent_index <- which.min(
  abs(
    autism_env$autism_data$baseprev *
      autism_env$autism_data$prev_factors - 0.01
  )
)
autism_one_percent_combined <- paired_combined_mutvar(
  autism_env$autism_data,
  autism_env$BurdenMLE_DN_models_autism,
  one_percent_index,
  "Combined Probands PTV",
  "Combined Probands Mis2",
  "Autism, 1% prevalence"
)
if (abs(autism_one_percent_combined$prevalence - 0.01) > 1e-10) {
  stop("The fitted autism models do not include prevalence 0.01.")
}
ddd_combined <- paired_combined_mutvar(
  ddd_env$kaplanis_data,
  ddd_env$BurdenMLE_DN_models_DDD,
  1L,
  "Proband PTV",
  "Proband Mis2",
  "Developmental disorders"
)
paired_combined <- rbind(
  autism_primary_combined,
  autism_one_percent_combined,
  ddd_combined
)

# The top-LOEUF fraction combines two fitted annotation strata. Calculate its
# interval by summing those strata within each existing gene-bootstrap
# replicate, preserving their covariance.
autism_model_index <- setNames(
  seq_along(autism_env$autism_data$loop_vars$names),
  autism_env$autism_data$loop_vars$names
)
ptv_index <- autism_model_index[["Combined Probands PTV"]]
ptv_model <- autism_env$BurdenMLE_DN_models_autism[[1]][[ptv_index]]
ptv_data <- get_genetic_data(
  ptv_index, autism_env$autism_data
)$genetic_data
loeuf1_indices <- which(
  colnames(ptv_model$features) %in% c("LOEUF1_mu1", "LOEUF1_mu2")
)
if (length(loeuf1_indices) != 2L) {
  stop("The PTV model does not contain both top-LOEUF annotation strata.")
}
loeuf1_fraction_boot <- vapply(
  seq_along(ptv_model$bootstrap_output$bootstrap_delta),
  function(iteration) {
    indices <- ptv_model$bootstrap_output$bootstrap_indices[, iteration]
    boot_model <- ptv_model
    boot_model$features <- ptv_model$features[indices, , drop = FALSE]
    boot_model$delta <- ptv_model$bootstrap_output$bootstrap_delta[[iteration]]
    annotation_mutvar <- annotation_mutvar_constants(
      boot_model, ptv_data[indices, , drop = FALSE]
    )
    sum(annotation_mutvar[loeuf1_indices]) / sum(annotation_mutvar)
  },
  numeric(1)
)
loeuf1_fraction_ci <- unname(quantile(
  loeuf1_fraction_boot, c(0.025, 0.975)
))

paired_combined_file <- file.path(
  output_dir, "paired_combined_mutational_variance.tsv"
)
write.table(
  paired_combined,
  paired_combined_file,
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Recompute input counts directly from the same fixed variant inputs used by
# set_up_asc.R. These are data audits, not fitted estimates.
variant_names <- c(
  "SPARK_iWES_v2_de_novo_calls_v1.1c_new_samples_only_2025-03-29.txt",
  "ASC_B17_B21_de_novos_v2.2_calls_2025-03-29.txt",
  "GeneDx_ASC_de_novos_GRCh38_v3.4_calls_2025-07-05.txt",
  "ASC_v17_published_autosomal_and_updated_XY_de_novo_calls_2025-03-29.txt",
  "ASC_B15_B16_published_autosomal_and_updated_XY_de_novo_calls_2025-03-29.txt",
  "SPARK_Pilot_GATK_published_autosomal_and_updated_XY_de_novo_calls_2025-03-29.txt",
  "SPARK_WES1_GATK_published_autosomal_and_updated_XY_de_novo_calls_2025-03-29.txt",
  "Kaplanis_DDD_ASD_de_novos_GRCh38_2025-04-09.txt"
)
variant_paths <- file.path(final_runs_dir, "inputs", "autism", "variants", variant_names)
if (any(!file.exists(variant_paths))) {
  stop("One or more fixed autism variant inputs are missing.")
}
variant_data <- lapply(
  variant_paths, read.delim, check.names = FALSE, stringsAsFactors = FALSE
)
variant_counts <- vapply(
  variant_data,
  function(data) {
    alpha_missense <- !is.na(data$am_pathogenicity) &
      data$am_pathogenicity >= 0.97
    mpc2 <- !is.na(data$MPC_v2) & data$MPC_v2 >= 2
    missense_score <- as.integer(alpha_missense) + as.integer(mpc2)
    proband <- data$Role == "Proband"
    c(
      PTV = sum((data$isPTV | data$isOS) & proband),
      PTV_indel = sum((data$isPTV | data$isOS) & data$isIndel & proband),
      Mis2 = sum(data$isMis & !data$isOS & missense_score == 2 & proband),
      Mis1 = sum(data$isMis & !data$isOS & missense_score == 1 & proband),
      Mis0 = sum(data$isMis & !data$isOS & missense_score == 0 & proband),
      Syn = sum(data$isSyn & !data$isOS & proband)
    )
  },
  numeric(6)
)
variant_counts <- rowSums(variant_counts)
variant_counts <- c(
  variant_counts,
  Missense = sum(variant_counts[c("Mis2", "Mis1", "Mis0")]),
  Total = sum(variant_counts[c("PTV", "Mis2", "Mis1", "Mis0", "Syn")])
)
manuscript_variant_counts <- c(
  Total = 46248,
  PTV = 6201,
  PTV_indel = 2944,
  Missense = 29045,
  Mis0 = 23585,
  Mis1 = 3315,
  Mis2 = 2145,
  Syn = 10992
)

cohort_n <- setNames(main_cohort$n_proband, main_cohort$dataset)
if (cohort_n[["All"]] != 38680) {
  stop("Main cohort N does not equal 38,680.")
}
subgroup_n <- function(cohort, group) {
  sum(cohort_subgroups$N_Proband[
    cohort_subgroups$Cohort == cohort &
      cohort_subgroups$Group %in% group
  ])
}
ddid_percent <- c(
  SPARK = 100 * subgroup_n(
    "SPARK", c("Male DDID", "Female DDID")
  ) / cohort_n[["SPARK"]],
  ASC = 100 * subgroup_n(
    "ASC", c("Male DDID", "Female DDID")
  ) / cohort_n[["ASC"]],
  GeneDx = 100 * subgroup_n(
    "GeneDx", c("Male DDID", "Female DDID")
  ) / cohort_n[["GeneDx"]]
)

constraint_length_file <- file.path(
  final_runs_dir, "inputs", "constraint",
  "gnomad.v2.1.1.lof_metrics.by_gene.txt"
)
ger_file <- file.path(final_runs_dir, "inputs", "reference", "GER_IDs.txt")
nc_file <- file.path(final_runs_dir, "inputs", "reference", "NC_IDs.txt")
constraint_length <- read.delim(
  constraint_length_file, check.names = FALSE, stringsAsFactors = FALSE
)
ger_genes <- read.delim(ger_file, header = FALSE, stringsAsFactors = FALSE)[[1]]
nc_genes <- read.delim(nc_file, header = FALSE, stringsAsFactors = FALSE)[[1]]
neither <- !constraint_length$gene_id %in% c(ger_genes, nc_genes)
median_cds <- c(
  Neither = median(constraint_length$cds_length[neither], na.rm = TRUE),
  GER = median(
    constraint_length$cds_length[
      constraint_length$gene_id %in% ger_genes
    ],
    na.rm = TRUE
  ),
  NC = median(
    constraint_length$cds_length[
      constraint_length$gene_id %in% nc_genes
    ],
    na.rm = TRUE
  )
)
cds_fold <- median_cds[c("GER", "NC")] / median_cds[["Neither"]]

overlap_pedigree_file <- file.path(
  final_runs_dir, "inputs", "autism", "pedigrees",
  paste0(
    "GeneDx_ASC_de_novos_GRCh38_v3.4_ped_HPO_terms_removed_",
    "2025-07-05_showing_Kaplanis_overlap.txt"
  )
)
overlap_pedigree <- read.delim(
  overlap_pedigree_file, check.names = FALSE, stringsAsFactors = FALSE
)
overlap_n <- length(unique(
  overlap_pedigree$Sample[overlap_pedigree$In_Kaplanis == 1]
))

ddd_model_index <- setNames(
  seq_along(ddd_env$kaplanis_data$loop_vars$names),
  ddd_env$kaplanis_data$loop_vars$names
)
ddd_excess_per_100 <- function(class) {
  index <- ddd_model_index[[paste("Proband", class)]]
  genetic_data <- get_genetic_data(index, ddd_env$kaplanis_data)$genetic_data
  scale <- if (class == "PTV") ddd_env$kaplanis_data$ptv_scale_factor else 1
  100 * scale *
    sum(genetic_data$case_count - genetic_data$expected_count) /
    genetic_data$N[1]
}
ddd_excess <- setNames(
  vapply(c("PTV", "Mis2", "Mis1", "Mis0"), ddd_excess_per_100, numeric(1)),
  c("PTV", "Mis2", "Mis1", "Mis0")
)

optimizer_current <- subset(optimizer_comparison, N == 38680)
if (!nrow(optimizer_current)) {
  stop("Optimizer comparison has no simulations at N = 38,680.")
}
optimizer_metrics <- c(
  total_runtime_speedup_median =
    median(optimizer_current$total_runtime_speedup),
  optimizer_runtime_speedup_median =
    median(optimizer_current$optimizer_runtime_speedup),
  absolute_loglik_difference_median =
    median(abs(optimizer_current$log_likelihood_difference)),
  absolute_loglik_difference_per_gene_max =
    max(abs(optimizer_current$log_likelihood_difference)) /
      length(unique(gene_table$gene_ID)),
  mutvar_estimate_correlation = cor(
    optimizer_current$est_mutvar_mixsqp,
    optimizer_current$est_mutvar_em
  ),
  mutvar_absolute_difference_median =
    median(abs(optimizer_current$estimated_mutvar_difference)),
  mixsqp_convergence_fraction =
    mean(optimizer_current$optimizer_converged_mixsqp)
)

format_percent_ci <- function(estimate, lower, upper, digits = 1) {
  sprintf(
    paste0("%.", digits, "f%% (95%% CI %.", digits, "f–%.", digits, "f%%)"),
    100 * estimate, 100 * lower, 100 * upper
  )
}
format_number_ci <- function(estimate, lower, upper, digits = 2) {
  sprintf(
    paste0("%.", digits, "f (95%% CI %.", digits, "f–%.", digits, "f)"),
    estimate, lower, upper
  )
}
format_count_ci <- function(estimate, lower, upper) {
  sprintf(
    "%d genes (95%% CI %d–%d)",
    round(estimate), round(lower), round(upper)
  )
}
format_plain <- function(value, digits = 2) {
  formatC(value, format = "f", digits = digits, big.mark = ",")
}

checklist <- data.frame()
add_row <- function(
    section, paragraph_anchor, estimand, stratum, estimate,
    lower = NA_real_, upper = NA_real_, unit = "",
    manuscript_ready, manuscript_snapshot = "",
    status = "ready", source_file, source_key = "", notes = "") {
  checklist <<- rbind(
    checklist,
    data.frame(
      order = nrow(checklist) + 1L,
      manuscript_section = section,
      paragraph_anchor = paragraph_anchor,
      estimand = estimand,
      stratum = stratum,
      estimate = as.numeric(estimate),
      ci_lower = as.numeric(lower),
      ci_upper = as.numeric(upper),
      unit = unit,
      manuscript_ready = manuscript_ready,
      manuscript_snapshot = manuscript_snapshot,
      update_status = status,
      source_file = relative_path(source_file),
      source_key = source_key,
      notes = notes,
      stringsAsFactors = FALSE
    )
  )
}

main_rr5 <- one_row(main_fraction, threshold == 5)
main_rr20 <- one_row(main_fraction, threshold == 20)
half_genes <- main_summary$half_mutvar_gene_count
scn2a <- one_row(gene_table, gene_name == "SCN2A")

# Abstract: repeated topline quantities remain separate rows so that the table
# follows manuscript order and can be checked paragraph by paragraph.
add_row(
  "Abstract", "Abstract", "Autism trio sample size", "All cohorts",
  cohort_n[["All"]], unit = "trios",
  manuscript_ready = "38,680 autism trios",
  manuscript_snapshot = "38,680 autism trios",
  source_file = main_cohort_file, source_key = "dataset=All"
)
add_row(
  "Abstract", "Abstract", "Combined PTV + Mis2 mutational variance",
  "Autism", autism_primary_combined$estimate,
  autism_primary_combined$lower, autism_primary_combined$upper,
  unit = "proportion",
  manuscript_ready = format_percent_ci(
    autism_primary_combined$estimate,
    autism_primary_combined$lower,
    autism_primary_combined$upper
  ),
  manuscript_snapshot = "3.2%",
  status = "update_value_and_add_CI",
  source_file = paired_combined_file,
  source_key = "Autism, primary prevalence"
)
add_row(
  "Abstract", "Abstract", "Fraction of cases with PTV/Mis2 RR > 5",
  "Autism", main_rr5$combined, main_rr5$combined_lower,
  main_rr5$combined_upper, unit = "proportion",
  manuscript_ready = format_percent_ci(
    main_rr5$combined, main_rr5$combined_lower, main_rr5$combined_upper
  ),
  manuscript_snapshot = "6.4%",
  status = "update_value_and_add_CI",
  source_file = main_fraction_file, source_key = "threshold=5"
)
add_row(
  "Abstract", "Abstract", "Genes explaining at least 50% of mutational variance",
  "Autism, PTV + Mis2", half_genes, unit = "genes",
  manuscript_ready = paste0(half_genes, " genes"),
  manuscript_snapshot = "14 genes",
  status = if (half_genes == 14) "matches" else "update_value",
  source_file = main_summary_file, source_key = "half_mutvar_gene_count"
)

# Results: data composition.
add_row(
  "Results", "Population and clinical genetic characterization",
  "Weighted autism prevalence", "All cohorts",
  main_summary$metadata$prevalence, unit = "proportion",
  manuscript_ready = paste0(
    format_plain(100 * main_summary$metadata$prevalence, 1), "%"
  ),
  manuscript_snapshot = "2.2%",
  source_file = main_summary_file, source_key = "metadata$prevalence"
)
for (cohort in c("All", "SPARK", "ASC", "GeneDx", "DDD")) {
  snapshot <- c(
    All = "38,680", SPARK = "18,060", ASC = "8,679",
    GeneDx = "10,747", DDD = "1,194"
  )[[cohort]]
  add_row(
    "Results", "Population and clinical genetic characterization",
    "Autism proband sample size", cohort, cohort_n[[cohort]], unit = "probands",
    manuscript_ready = paste0("N = ", format(cohort_n[[cohort]], big.mark = ",")),
    manuscript_snapshot = paste0("N = ", snapshot),
    source_file = main_cohort_file, source_key = paste0("dataset=", cohort)
  )
}
for (entry in c("Total", "PTV", "PTV_indel", "Missense", "Mis0", "Mis1", "Mis2", "Syn")) {
  label <- c(
    Total = "All coding de novo variants",
    PTV = "Protein-truncating variants",
    PTV_indel = "PTV indels",
    Missense = "Missense variants",
    Mis0 = "Mis0 variants",
    Mis1 = "Mis1 variants",
    Mis2 = "Mis2 variants",
    Syn = "Synonymous variants"
  )[[entry]]
  add_row(
    "Results", "Population and clinical genetic characterization",
    label, "Autism probands", variant_counts[[entry]], unit = "variants",
    manuscript_ready = format(variant_counts[[entry]], big.mark = ","),
    manuscript_snapshot = format(
      manuscript_variant_counts[[entry]], big.mark = ","
    ),
    status = if (
      variant_counts[[entry]] == manuscript_variant_counts[[entry]]
    ) "matches" else "update_value",
    source_file = variant_paths[1],
    source_key = paste0("fixed input set; ", entry),
    notes = "Recomputed across all eight fixed variant input files."
  )
}

# Mutational-variance estimates and intervals for each class.
for (class in c("PTV", "Mis2", "Mis1", "Mis0", "Syn")) {
  row <- one_row(main_mutvar, variant_class == class)
  digits <- if (class %in% c("Mis0", "Syn")) 2 else 1
  add_row(
    "Results", "Population and clinical genetic characterization",
    "Mutational variance", class, row$estimate, row$lower, row$upper,
    unit = "proportion",
    manuscript_ready = format_percent_ci(
      row$estimate, row$lower, row$upper, digits
    ),
    manuscript_snapshot = c(
      PTV = "2.7%", Mis2 = "0.6%", Mis1 = "0.1%",
      Mis0 = "0.06%", Syn = "0.04%"
    )[[class]],
    status = "update_value_and_add_CI",
    source_file = main_mutvar_file,
    source_key = paste0("variant_class=", class)
  )
}
add_row(
  "Results", "Population and clinical genetic characterization",
  "Combined PTV + Mis2 mutational variance", "Autism",
  autism_primary_combined$estimate,
  autism_primary_combined$lower,
  autism_primary_combined$upper,
  unit = "proportion",
  manuscript_ready = format_percent_ci(
    autism_primary_combined$estimate,
    autism_primary_combined$lower,
    autism_primary_combined$upper
  ),
  manuscript_snapshot = "3.2%",
  status = "update_value_and_add_CI",
  source_file = paired_combined_file,
  source_key = "Autism, primary prevalence",
  notes = "Paired gene-bootstrap interval."
)

paternal_spark <- one_row(
  paternal_summary, group == "SPARK autistic probands"
)
paternal_us <- one_row(paternal_summary, group == "U.S. births (2024)")
add_row(
  "Results", "Population and clinical genetic characterization",
  "Mean paternal-age difference", "SPARK minus U.S.",
  paternal_test$mean_difference_years,
  paternal_test$ci_lower, paternal_test$ci_upper,
  unit = "years",
  manuscript_ready = format_number_ci(
    paternal_test$mean_difference_years,
    paternal_test$ci_lower,
    paternal_test$ci_upper,
    2
  ),
  manuscript_snapshot = "~5 months",
  source_file = paternal_test_file,
  source_key = "Welch mean comparison",
  notes = paste0(
    "SPARK mean ", format_plain(paternal_spark$mean, 2),
    "; U.S. mean ", format_plain(paternal_us$mean, 2), "."
  )
)
add_row(
  "Results", "Population and clinical genetic characterization",
  "Paternal-age standardized mean difference", "SPARK minus U.S.",
  paternal_test$cohens_d, unit = "Cohen's d",
  manuscript_ready = paste0("Cohen’s d = ", format_plain(paternal_test$cohens_d, 2)),
  manuscript_snapshot = "Cohen’s d = 0.06",
  source_file = paternal_test_file, source_key = "cohens_d"
)
paternal_combined <- one_row(
  paternal_bias, variant_class == "Combined PTV + Mis2"
)
add_row(
  "Results", "Population and clinical genetic characterization",
  "Paternal-age-associated mutational-variance change",
  "Combined PTV + Mis2",
  -paternal_combined$absolute_change,
  unit = "absolute mutational variance",
  manuscript_ready = paste0(
    "estimated upward bias of ",
    formatC(-paternal_combined$absolute_change, format = "f", digits = 5)
  ),
  manuscript_snapshot = "upward bias of 0.00019",
  source_file = paternal_bias_file,
  source_key = "Combined PTV + Mis2",
  notes = "Magnitude is reported as upward bias; the sensitivity table stores the downward correction."
)
add_row(
  "Results", "Population and clinical genetic characterization",
  "Combined PTV + Mis2 mutational variance", "Autism; prevalence 1%",
  autism_one_percent_combined$estimate,
  autism_one_percent_combined$lower,
  autism_one_percent_combined$upper,
  unit = "proportion",
  manuscript_ready = format_percent_ci(
    autism_one_percent_combined$estimate,
    autism_one_percent_combined$lower,
    autism_one_percent_combined$upper
  ),
  manuscript_snapshot = "1.4%",
  status = "add_CI",
  source_file = paired_combined_file,
  source_key = "Autism, 1% prevalence",
  notes = "Paired gene-bootstrap interval."
)
add_row(
  "Results", "Population and clinical genetic characterization",
  "Spearman correlation with TADA Bayes factor",
  "Posterior-mean PTV rate ratio",
  tada$estimate[grepl("posterior mean", tada$statistic)],
  unit = "Spearman rho",
  manuscript_ready = paste0(
    "Spearman’s ρ = ",
    format_plain(tada$estimate[grepl("posterior mean", tada$statistic)], 2)
  ),
  manuscript_snapshot = "positively correlated",
  status = "optional_add_numeric_result",
  source_file = tada_file,
  source_key = "posterior mean vs BF_PTV"
)
add_row(
  "Results", "Population and clinical genetic characterization",
  "Korean WGS affected-offspring sample size", "Kim et al.",
  unique(korean_results$n_probands), unit = "probands",
  manuscript_ready = paste0(
    format(unique(korean_results$n_probands), big.mark = ","),
    " affected offspring"
  ),
  manuscript_snapshot = "680 affected offspring",
  source_file = korean_results_file, source_key = "n_probands"
)
for (class in c("PTV", "Pooled missense")) {
  row <- one_row(korean_results, variant_class == class)
  digits <- if (class == "PTV") 4 else 4
  add_row(
    "Results", "Population and clinical genetic characterization",
    "Korean WGS mutational variance", class,
    row$mutational_variance, row$mutvar_ci_lower, row$mutvar_ci_upper,
    unit = "proportion",
    manuscript_ready = format_number_ci(
      row$mutational_variance, row$mutvar_ci_lower, row$mutvar_ci_upper, digits
    ),
    manuscript_snapshot = if (class == "PTV") {
      "0.0009 (95% CI 0.0004–0.0443)"
    } else {
      "0.0215 (95% CI 0.0046–0.0559)"
    },
    source_file = korean_results_file,
    source_key = paste0("variant_class=", class)
  )
}

# Polygenicity and annotations.
add_row(
  "Results", "Polygenicity", "Genes explaining at least 50% of mutational variance",
  "Autism, PTV + Mis2", half_genes, unit = "genes",
  manuscript_ready = paste0(half_genes, " genes"),
  manuscript_snapshot = "14 genes",
  status = if (half_genes == 14) "matches" else "update_value",
  source_file = main_summary_file, source_key = "half_mutvar_gene_count"
)
near_all_genes <- main_summary$polygenicity$number_genes[
  which(main_summary$polygenicity$cumulative_mutvar >= 0.95)[1]
]
add_row(
  "Results", "Polygenicity", "Genes explaining at least 95% of mutational variance",
  "Autism, PTV + Mis2", near_all_genes, unit = "genes",
  manuscript_ready = paste0(
    "approximately ", format(round(near_all_genes, -2), big.mark = ","),
    " genes"
  ),
  manuscript_snapshot = "approximately 1000 genes",
  status = "update_value",
  source_file = main_summary_file, source_key = "polygenicity cumulative_mutvar>=0.95",
  notes = "A 95% threshold operationalizes the manuscript phrase 'nearly 100%'."
)
add_row(
  "Results", "Polygenicity", "Fraction of mutational variance explained",
  "SCN2A", scn2a$gene_FractionMutVar, unit = "proportion",
  manuscript_ready = paste0(
    format_plain(100 * scn2a$gene_FractionMutVar, 1), "%"
  ),
  manuscript_snapshot = "8.7%",
  status = "update_value",
  source_file = gene_table_file, source_key = "gene_name=SCN2A"
)
top_constraint_fraction <- sum(
  main_summary$enrichment$fraction_mutvar[
    main_summary$enrichment$annotation %in% c("LOEUF1_mu1", "LOEUF1_mu2")
  ]
)
add_row(
  "Results", "Annotation enrichment",
  "Fraction of PTV mutational variance in top LOEUF quintile",
  "LOEUF quintile 1", top_constraint_fraction,
  loeuf1_fraction_ci[1], loeuf1_fraction_ci[2],
  unit = "proportion",
  manuscript_ready = format_percent_ci(
    top_constraint_fraction,
    loeuf1_fraction_ci[1],
    loeuf1_fraction_ci[2]
  ),
  manuscript_snapshot = "89%",
  status = "update_value_and_add_CI",
  source_file = autism_model_file,
  source_key = "Combined Probands PTV; LOEUF1_mu1 + LOEUF1_mu2",
  notes = "Percentile interval from the existing gene-bootstrap replicates."
)
for (gene_set in c("GER", "NC")) {
  add_row(
    "Results", "Annotation enrichment",
    "Median coding-sequence-length ratio",
    if (gene_set == "GER") {
      "Gene-expression regulation vs neither"
    } else {
      "Neuronal communication vs neither"
    },
    cds_fold[[gene_set]], unit = "fold",
    manuscript_ready = paste0(format_plain(cds_fold[[gene_set]], 1), "x"),
    manuscript_snapshot = if (gene_set == "GER") "1.2x" else "1.3x",
    status = "check_updated_value",
    source_file = constraint_length_file,
    source_key = gene_set,
    notes = "Uses the same gene sets and median CDS calculation as Supplementary Figure 8."
  )
}

# Fraction of cases, effective penetrance, and effective rate ratio.
for (rr_threshold in c(5, 20)) {
  row <- one_row(main_fraction, threshold == rr_threshold)
  add_row(
    "Results", "Fraction of cases and effective penetrance",
    "Fraction of cases with PTV/Mis2 above rate-ratio threshold",
    paste0("RR > ", rr_threshold),
    row$combined, row$combined_lower, row$combined_upper,
    unit = "proportion",
    manuscript_ready = format_percent_ci(
      row$combined, row$combined_lower, row$combined_upper
    ),
    manuscript_snapshot = if (rr_threshold == 5) "6.4%" else "2.6%",
    status = "update_value_and_add_CI",
    source_file = main_fraction_file,
    source_key = paste0("threshold=", rr_threshold)
  )
}
for (rr_threshold in c(5, 20)) {
  threshold_penetrance <- main_summary$metadata$prevalence * rr_threshold
  add_row(
    "Results", "Fraction of cases and effective penetrance",
    "Diagnosis probability corresponding to rate-ratio threshold",
    paste0("RR = ", rr_threshold),
    threshold_penetrance, unit = "probability",
    manuscript_ready = format_plain(threshold_penetrance, 2),
    manuscript_snapshot = if (rr_threshold == 5) "0.11" else "0.44",
    source_file = main_summary_file,
    source_key = paste0("prevalence * ", rr_threshold)
  )
}
for (class in c("PTV", "Mis2", "Mis1")) {
  row <- one_row(main_penetrance, variant_class == class)
  add_row(
    "Results", "Fraction of cases and effective penetrance",
    "Effective penetrance", class,
    row$estimate, row$lower, row$upper, unit = "probability",
    manuscript_ready = format_number_ci(
      row$estimate, row$lower, row$upper, 2
    ),
    manuscript_snapshot = c(PTV = "0.48", Mis2 = "0.25", Mis1 = "0.10")[[class]],
    status = "update_value_and_add_CI",
    source_file = main_penetrance_file,
    source_key = paste0("variant_class=", class)
  )
}
for (class in c("PTV", "Mis2", "Mis1")) {
  row <- one_row(
    autism_dd_rr,
    diagnosis == "Autism" & variant_class == class
  )
  add_row(
    "Results", "Fraction of cases and effective penetrance",
    "Effective rate ratio", class,
    row$estimate, row$lower, row$upper, unit = "rate ratio",
    manuscript_ready = format_number_ci(
      row$estimate, row$lower, row$upper, 1
    ),
    manuscript_snapshot = c(PTV = "22", Mis2 = "11", Mis1 = "4.7")[[class]],
    status = "update_value_and_add_CI",
    source_file = autism_dd_rr_file,
    source_key = paste0("Autism; variant_class=", class)
  )
}
for (cohort in c("GeneDx", "ASC", "SPARK")) {
  add_row(
    "Results", "Cohort, sex, and DDID heterogeneity",
    "Proportion of probands with DDID", cohort,
    ddid_percent[[cohort]] / 100, unit = "proportion",
    manuscript_ready = paste0(format_plain(ddid_percent[[cohort]], 1), "%"),
    manuscript_snapshot = c(
      GeneDx = "42%", ASC = "23%", SPARK = "16%"
    )[[cohort]],
    status = if (
      cohort != "GeneDx" &&
        abs(ddid_percent[[cohort]] - c(ASC = 23, SPARK = 16)[[cohort]]) < 1
    ) "matches_rounded" else "update_value",
    source_file = cohort_subgroup_file,
    source_key = paste0(cohort, "; Male DDID + Female DDID"),
    notes = "Computed from the current cohort subgroup Ns."
  )
}

# Future gene-discovery estimates and intervals.
for (rr_threshold in c(2, 5, 10, 20)) {
  row <- one_row(
    autism_dd_genes,
    diagnosis == "Autism" & threshold == rr_threshold
  )
  add_row(
    "Results", "Future autism gene discovery",
    "Number of PTV-associated autism genes",
    paste0("RR > ", rr_threshold),
    row$estimate, row$lower, row$upper, unit = "genes",
    manuscript_ready = format_count_ci(row$estimate, row$lower, row$upper),
    manuscript_snapshot = c(
      `2` = "526", `5` = "263", `10` = "125", `20` = "47"
    )[[as.character(rr_threshold)]],
    status = "update_value_and_add_CI",
    source_file = autism_dd_genes_file,
    source_key = paste0("Autism; threshold=", rr_threshold),
    notes = "Interval for the estimated total number of genes."
  )
}
gene_probability <- main_summary$ptv_gene_rr_probability
probability_columns <- paste0("rr_greater_", c(2, 5, 10, 20))
for (fdr_threshold in c(0.001, 0.05)) {
  fdr_label <- paste0("ASC 2026 FDR < ", fdr_threshold)
  significant_genes <- posterior_summary$tada$Gene_ID[
    !posterior_summary$tada$Flag &
      !is.na(posterior_summary$tada$FDR) &
      posterior_summary$tada$FDR < fdr_threshold
  ]
  matched <- gene_probability$gene_id %in% significant_genes
  discovered_counts <- colSums(
    gene_probability[matched, probability_columns, drop = FALSE]
  )
  for (rr_threshold in c(2, 5, 10, 20)) {
    count <- discovered_counts[[paste0("rr_greater_", rr_threshold)]]
    add_row(
      "Results", "Future autism gene discovery",
      "Previously discovered genes",
      paste0(fdr_label, "; RR > ", rr_threshold),
      count, unit = "genes",
      manuscript_ready = paste0(round(count), " genes"),
      manuscript_snapshot = "",
      status = "update_from_ASC_2026",
      source_file = posterior_summary_file,
      source_key = paste0(fdr_label, "; RR > ", rr_threshold),
      notes = "Posterior expected count within the current ASC-significant gene set."
    )
  }
}

# Developmental-disorder comparison.
ddd_n <- ddd_env$kaplanis_data$loop_vars$N_subset[[
  ddd_model_index[["Proband PTV"]]
]]
add_row(
  "Results", "Comparison with severe developmental disorders",
  "Developmental-disorder trio sample size", "Kaplanis et al.",
  ddd_n, unit = "trios",
  manuscript_ready = paste0(format(ddd_n, big.mark = ","), " trios"),
  manuscript_snapshot = "31,058 trios",
  source_file = ddd_model_file, source_key = "Proband PTV N"
)
for (class in c("PTV", "Mis2", "Mis1", "Mis0", "Syn")) {
  row <- one_row(
    autism_dd_mutvar,
    diagnosis == "DD" & variant_class == class
  )
  digits <- if (class %in% c("Mis0", "Syn")) 2 else 1
  add_row(
    "Results", "Comparison with severe developmental disorders",
    "Developmental-disorder mutational variance", class,
    row$estimate, row$lower, row$upper, unit = "proportion",
    manuscript_ready = format_percent_ci(
      row$estimate, row$lower, row$upper, digits
    ),
    manuscript_snapshot = c(
      PTV = "6.7%", Mis2 = "1.7%", Mis1 = "0.5%",
      Mis0 = "0.1%", Syn = "0.02%"
    )[[class]],
    status = "update_value_and_add_CI",
    source_file = autism_dd_mutvar_file,
    source_key = paste0("DD; variant_class=", class)
  )
}
add_row(
  "Results", "Comparison with severe developmental disorders",
  "Combined PTV + Mis2 mutational variance", "Developmental disorders",
  ddd_combined$estimate, ddd_combined$lower, ddd_combined$upper,
  unit = "proportion",
  manuscript_ready = format_percent_ci(
    ddd_combined$estimate, ddd_combined$lower, ddd_combined$upper
  ),
  manuscript_snapshot = "8.4%",
  status = "update_value_and_add_CI",
  source_file = paired_combined_file,
  source_key = "Developmental disorders",
  notes = "Paired gene-bootstrap interval."
)
add_row(
  "Results", "Comparison with severe developmental disorders",
  "DD-to-autism combined mutational-variance ratio", "PTV + Mis2",
  ddd_combined$estimate / autism_primary_combined$estimate,
  unit = "fold",
  manuscript_ready = paste0(
    format_plain(ddd_combined$estimate / autism_primary_combined$estimate, 1),
    "-fold"
  ),
  manuscript_snapshot = "approximately 2.5x",
  status = "update_value",
  source_file = paired_combined_file,
  source_key = "DD / Autism primary"
)
ddd_rr5 <- one_row(
  autism_dd_fraction,
  diagnosis == "DD" & threshold == 5
)
add_row(
  "Results", "Comparison with severe developmental disorders",
  "Fraction of DD cases with PTV/Mis2 RR > 5", "Developmental disorders",
  ddd_rr5$estimate, ddd_rr5$lower, ddd_rr5$upper,
  unit = "proportion",
  manuscript_ready = format_percent_ci(
    ddd_rr5$estimate, ddd_rr5$lower, ddd_rr5$upper
  ),
  manuscript_snapshot = "16%",
  status = "update_value_and_add_CI",
  source_file = autism_dd_fraction_file,
  source_key = "DD; threshold=5"
)
ddd_ptv_rr <- one_row(
  autism_dd_rr,
  diagnosis == "DD" & variant_class == "PTV"
)
add_row(
  "Results", "Comparison with severe developmental disorders",
  "Effective rate ratio", "DD PTV",
  ddd_ptv_rr$estimate, ddd_ptv_rr$lower, ddd_ptv_rr$upper,
  unit = "rate ratio",
  manuscript_ready = format_number_ci(
    ddd_ptv_rr$estimate, ddd_ptv_rr$lower, ddd_ptv_rr$upper, 1
  ),
  manuscript_snapshot = "64",
  status = "add_CI",
  source_file = autism_dd_rr_file,
  source_key = "DD; variant_class=PTV"
)
add_row(
  "Results", "Comparison with severe developmental disorders",
  "Autism/DD overlapping samples", "GeneDx/Kaplanis overlap",
  overlap_n, unit = "samples",
  manuscript_ready = paste0(format(overlap_n, big.mark = ","), " samples"),
  manuscript_snapshot = "3,543 samples",
  status = if (overlap_n == 3543) "matches" else "update_value",
  source_file = overlap_pedigree_file,
  source_key = "unique Sample with In_Kaplanis=1"
)
add_row(
  "Results", "Comparison with severe developmental disorders",
  "Excess variants per 100 cases", "DD PTV",
  ddd_excess[["PTV"]], unit = "variants per 100 cases",
  manuscript_ready = format_plain(ddd_excess[["PTV"]], 1),
  manuscript_snapshot = "10.0",
  status = "update_value",
  source_file = ddd_model_file, source_key = "Proband PTV observed-expected"
)
add_row(
  "Results", "Comparison with severe developmental disorders",
  "Excess variants per 100 cases", "DD missense (Mis0+Mis1+Mis2)",
  sum(ddd_excess[c("Mis0", "Mis1", "Mis2")]),
  unit = "variants per 100 cases",
  manuscript_ready = format_plain(
    sum(ddd_excess[c("Mis0", "Mis1", "Mis2")]), 1
  ),
  manuscript_snapshot = "25.9",
  status = "update_value",
  source_file = ddd_model_file,
  source_key = "Proband Mis0+Mis1+Mis2 observed-expected"
)

# Discussion repetitions are intentionally retained as independent checklist
# rows because these values must also be updated in prose.
for (item in list(
  list("Combined PTV + Mis2 mutational variance", autism_primary_combined,
       "3.2%", format_percent_ci(
         autism_primary_combined$estimate,
         autism_primary_combined$lower,
         autism_primary_combined$upper
       )),
  list("Fraction of cases with PTV/Mis2 RR > 5",
       data.frame(
         estimate = main_rr5$combined,
         lower = main_rr5$combined_lower,
         upper = main_rr5$combined_upper
       ),
       "6.4%", format_percent_ci(
         main_rr5$combined, main_rr5$combined_lower, main_rr5$combined_upper
       ))
)) {
  add_row(
    "Discussion", "Discussion summary",
    item[[1]], "Autism",
    item[[2]]$estimate, item[[2]]$lower, item[[2]]$upper,
    unit = "proportion",
    manuscript_ready = item[[4]],
    manuscript_snapshot = item[[3]],
    status = "update_repeated_value",
    source_file = if (grepl("Fraction", item[[1]])) {
      main_fraction_file
    } else {
      paired_combined_file
    },
    source_key = "Repeated topline estimate"
  )
}
add_row(
  "Discussion", "Discussion summary",
  "Genes explaining at least 50% of mutational variance",
  "Autism, PTV + Mis2", half_genes, unit = "genes",
  manuscript_ready = paste0(half_genes, " genes"),
  manuscript_snapshot = "14 genes",
  status = if (half_genes == 14) "matches" else "update_repeated_value",
  source_file = main_summary_file, source_key = "half_mutvar_gene_count"
)

# Proposed new Methods diagnostics. These are deliberately based on paired
# simulations at the observed study N rather than on unmatched legacy models.
optimizer_rows <- list(
  list(
    "Median end-to-end fitting speedup", optimizer_metrics[[
      "total_runtime_speedup_median"
    ]], "fold", paste0(
      format_plain(optimizer_metrics[["total_runtime_speedup_median"]], 2),
      "-fold faster"
    )
  ),
  list(
    "Median optimizer-only speedup", optimizer_metrics[[
      "optimizer_runtime_speedup_median"
    ]], "fold", paste0(
      format_plain(optimizer_metrics[["optimizer_runtime_speedup_median"]], 1),
      "-fold faster"
    )
  ),
  list(
    "Median absolute achieved log-likelihood difference",
    optimizer_metrics[["absolute_loglik_difference_median"]],
    "log-likelihood units",
    format_plain(
      optimizer_metrics[["absolute_loglik_difference_median"]], 2
    )
  ),
  list(
    "Maximum absolute achieved log-likelihood difference per gene",
    optimizer_metrics[["absolute_loglik_difference_per_gene_max"]],
    "log-likelihood units per gene",
    formatC(
      optimizer_metrics[["absolute_loglik_difference_per_gene_max"]],
      format = "f", digits = 6
    )
  ),
  list(
    "Correlation of MixSQP and EM mutational-variance estimates",
    optimizer_metrics[["mutvar_estimate_correlation"]],
    "Pearson correlation",
    format_plain(
      optimizer_metrics[["mutvar_estimate_correlation"]], 5
    )
  ),
  list(
    "Median absolute MixSQP-EM mutational-variance difference",
    optimizer_metrics[["mutvar_absolute_difference_median"]],
    "mutational variance",
    formatC(
      optimizer_metrics[["mutvar_absolute_difference_median"]],
      format = "f", digits = 6
    )
  ),
  list(
    "MixSQP convergence fraction", optimizer_metrics[[
      "mixsqp_convergence_fraction"
    ]], "proportion", paste0(
      format_plain(
        100 * optimizer_metrics[["mixsqp_convergence_fraction"]], 1
      ),
      "%"
    )
  )
)
for (item in optimizer_rows) {
  add_row(
    "Methods", "Proposed MixSQP implementation paragraph",
    item[[1]], "Paired simulations; N = 38,680",
    item[[2]], unit = item[[3]],
    manuscript_ready = item[[4]],
    manuscript_snapshot = "",
    status = "new_methods_text",
    source_file = simulation_diagnostic_file,
    source_key = "N=38680",
    notes = paste0(
      "Candidate topline diagnostic for text accompanying ",
      "Supplementary Figure 17; 400 paired simulations."
    )
  )
}

# Basic analysis specifications that are easy to lose during manuscript edits.
for (specification in list(
  list("Number of uniform mixture components", 10, "components", "10"),
  list("Maximum fitted rate ratio", 100, "rate ratio", "100"),
  list("Poisson-Uniform evaluation grid size", 10, "grid points", "10"),
  list("Gene-bootstrap iterations", main_summary$metadata$n_boot,
       "iterations", as.character(main_summary$metadata$n_boot)),
  list("gnomAD calibration mixture components", 31,
       "components", "31")
)) {
  add_row(
    "Methods", "Implementation details",
    specification[[1]], "Primary analysis",
    specification[[2]], unit = specification[[3]],
    manuscript_ready = specification[[4]],
    manuscript_snapshot = specification[[4]],
    source_file = if (
      specification[[1]] == "Gene-bootstrap iterations"
    ) main_summary_file else file.path(repo_dir, "analysis", "scripts", "set_up_asc.R"),
    source_key = "Analysis specification"
  )
}

incomplete_ci <- xor(
  is.finite(checklist$ci_lower),
  is.finite(checklist$ci_upper)
)
if (any(incomplete_ci)) {
  bad <- checklist$estimand[incomplete_ci]
  stop(
    "Rows have only one confidence-interval bound: ",
    paste(unique(bad), collapse = ", ")
  )
}
complete_ci <- is.finite(checklist$ci_lower) & is.finite(checklist$ci_upper)
invalid_ci <- complete_ci & (
  checklist$ci_lower > checklist$estimate |
    checklist$estimate > checklist$ci_upper
)
if (any(invalid_ci)) {
  stop(
    "Point estimates fall outside their confidence intervals: ",
    paste(unique(checklist$estimand[invalid_ci]), collapse = ", ")
  )
}
required_gene_thresholds <- paste0("RR > ", c(2, 5, 10, 20))
reported_gene_thresholds <- checklist$stratum[
  checklist$manuscript_section == "Results" &
    checklist$estimand == "Number of PTV-associated autism genes"
]
if (!all(required_gene_thresholds %in% reported_gene_thresholds)) {
  stop("The future-gene checklist is incomplete.")
}

checklist_file <- file.path(output_dir, "ManuscriptEstimateChecklist.tsv")
write.table(
  checklist, checklist_file,
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

supplementary_inventory <- data.frame(
  figure = paste0("Supplementary Figure ", 1:16),
  file = c(
    "SupplementaryFigure1_FullSimulationResults.pdf",
    "SupplementaryFigure2_SiblingNegativeControl.pdf",
    "SupplementaryFigure3_CESGeneExclusion.pdf",
    "SupplementaryFigure4_PaternalAge.pdf",
    "SupplementaryFigure5_PrevalenceSensitivity.pdf",
    "SupplementaryFigure6_TADAComparison.pdf",
    "SupplementaryFigure7_KoreanWGS.pdf",
    "SupplementaryFigure8_GeneSetCodingLength.pdf",
    "SupplementaryFigure9_CodingLengthByRateRatio.pdf",
    "SupplementaryFigure10_CohortHeterogeneity.pdf",
    "SupplementaryFigure11_DDIDComparison.pdf",
    "SupplementaryFigure12_DDIDMaximallyStratified.pdf",
    "SupplementaryFigure13_DDIDForecasting.pdf",
    "SupplementaryFigure14_DDPrevalenceSensitivity.pdf",
    "SupplementaryFigure15_AutismDDNoOverlap.pdf",
    "SupplementaryFigure16_OptimizerComparison.pdf"
  ),
  manuscript_location = c(
    "Simulation Results",
    "Autism mutational variance Results",
    "Autism mutational variance Results",
    "Autism mutational variance Results",
    "Autism mutational variance Results",
    "Autism mutational variance Results",
    "Autism mutational variance Results",
    "Annotation Results",
    "Annotation Results",
    "Cohort heterogeneity Results",
    "DDID Results",
    "DDID Results",
    "Future gene discovery Results",
    "Developmental-disorder Results",
    "Developmental-disorder Results",
    "MixSQP implementation Methods"
  )
)
write.table(
  supplementary_inventory,
  file.path(output_dir, "SupplementaryFigureInventory.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("Manuscript estimate checklist written to:\n")
cat(" ", checklist_file, "\n")
cat("Rows:", nrow(checklist), "\n")
cat(
  "Rows flagged for a numerical/text update:",
  sum(grepl("update|new_methods", checklist$update_status)),
  "\n"
)
