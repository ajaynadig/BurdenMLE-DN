# Build manuscript supplementary tables from the current MixSQP model fits and
# derived summaries.

suppressPackageStartupMessages({
  library(data.table)
  library(SummarizedExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit)) return(sub(paste0("^", flag, "="), "", hit[1]))
  i <- match(flag, args)
  if (!is.na(i)) return(args[i + 1])
  default
}

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo_dir <- dirname(dirname(dirname(normalizePath(script_file, mustWork = TRUE))))
final_runs_dir <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT", unset = file.path(repo_dir, "analysis")), mustWork = FALSE)
derived_dir <- file.path(final_runs_dir, "outputs", "derived")
output_dir <- get_arg(
  "--output-dir", file.path(final_runs_dir, "outputs", "tables", "supplementary")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

model_manifest <- get_arg("--model-manifest", NA_character_)
legacy_autism_file <- get_arg("--legacy-autism-model-file", NA_character_)
legacy_ddd_file <- get_arg("--legacy-ddd-model-file", NA_character_)

for (source_file in c(
  "BurdenMLE_DN.R", "estimate_mutvar.R", "io.R", "likelihoods.R",
  "model.R", "secondary_analysis_functions.R"
)) source(if (identical(source_file, "secondary_analysis_functions.R")) file.path(repo_dir, "analysis", "scripts", source_file) else file.path(repo_dir, "R", source_file))
source(file.path(repo_dir, "analysis", "scripts", "model_artifacts.R"))
autism_info <- load_model_artifact(
  model_manifest, "autism", envir = environment(),
  legacy_path = legacy_autism_file
)
ddd_info <- load_model_artifact(
  model_manifest, "ddd", envir = environment(),
  legacy_path = legacy_ddd_file
)
autism_model_file <- autism_info$path
ddd_model_file <- ddd_info$path

write_tsv <- function(x, filename) {
  write.table(x, file.path(output_dir, filename), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
}

variant_class <- function(name) {
  classes <- c("PTV", "Mis2", "Mis1", "Mis0", "Syn")
  hit <- classes[vapply(classes, function(x) grepl(x, name, fixed = TRUE), logical(1))]
  if (length(hit)) hit[1] else NA_character_
}

model_summary_table <- function(data, fits, keep_names) {
  annotations <- c("LOEUF1_mu1", "LOEUF1_mu2", paste0("LOEUF", 2:5))
  rows <- list()
  row_number <- 0L
  for (prevalence_index in seq_along(data$prev_factors)) {
    for (model_index in which(data$loop_vars$names %in% keep_names)) {
      model <- fits[[prevalence_index]][[model_index]]
      if (length(model) == 1L && is.na(model)) next
      row_number <- row_number + 1L
      name <- data$loop_vars$names[model_index]
      class <- variant_class(name)
      ptv_scale <- if (identical(class, "PTV")) data$ptv_scale_factor else 1
      row <- data.frame(
        `Model Name` = name,
        `Variant Class` = class,
        Prevalence = data$loop_vars$prevalences[model_index] *
          data$prev_factors[prevalence_index],
        Role = if (grepl("Sibling", name)) "Sibling" else "Proband",
        MutVar = model$mutvar_output$total_mutvar * ptv_scale,
        MutVar_lower95CI = model$mutvar_output$mutvar_CI[1] * ptv_scale,
        MutVar_upper95CI = model$mutvar_output$mutvar_CI[2] * ptv_scale,
        EffectivePenetrance = model$penetrance$effective_penetrance,
        EffectivePenetrance_Lower95CI = model$penetrance$effective_penetrance_CI[1],
        EffectivePenetrance_Upper95CI = model$penetrance$effective_penetrance_CI[2],
        check.names = FALSE
      )
      for (annotation_index in seq_along(annotations)) {
        annotation <- annotations[annotation_index]
        row[[paste0("Enrichment_", annotation)]] <-
          model$mutvar_output$enrichment[annotation_index]
        row[[paste0("Enrichment_Lower95CI_", annotation)]] <-
          model$mutvar_output$enrich_CI[1, annotation_index]
        row[[paste0("Enrichment_Upper95CI_", annotation)]] <-
          model$mutvar_output$enrich_CI[2, annotation_index]
        row[[paste0("FractionMutVar_", annotation)]] <-
          model$mutvar_output$frac_mutvar[annotation_index]
        row[[paste0("FractionExpected_", annotation)]] <-
          model$mutvar_output$frac_expected[annotation_index]
      }
      rows[[row_number]] <- row
    }
  }
  do.call(rbind, rows)
}

# Table 1: simulation estimates for the optimizer used in the manuscript run.
simulation_table_file <- file.path(
  output_dir, "SupplementaryTable1_SimulationResults.tsv"
)
if (!file.exists(simulation_table_file)) stop("Current MixSQP simulation table is missing.")
table1 <- fread(simulation_table_file, data.table = FALSE)
write_tsv(table1, "SupplementaryTable1_SimulationResults.tsv")

# Table 2: the constrained-PTV count table requested during revision. This is
# the expanded counterpart to the cohort summary in Figure 2A.
groups <- c(
  "All", "Male DDID", "Male Non-DDID",
  "Female DDID", "Female Non-DDID"
)
cohort_model_names <- list(
  All = c(
    "Combined Probands PTV", "Male DDID PTV", "Male Non-DDID PTV",
    "Female DDID PTV", "Female Non-DDID PTV"
  ),
  SPARK = paste(
    c("Combined", "Male DDID", "Male Non-DDID",
      "Female DDID", "Female Non-DDID"),
    "SPARK PTV"
  ),
  ASC = paste(
    c("Combined", "Male DDID", "Male Non-DDID",
      "Female DDID", "Female Non-DDID"),
    "ASC PTV"
  ),
  GeneDx = paste(
    c("Combined", "Male DDID", "Male Non-DDID",
      "Female DDID", "Female Non-DDID"),
    "GeneDx PTV"
  )
)
constrained_ptv_row <- function(cohort, group, model_data) {
  constrained <- model_data$features[, 1] == 1 |
    model_data$features[, 2] == 1
  observed <- sum(model_data$genetic_data$case_count[constrained])
  expected <- sum(model_data$genetic_data$expected_count[constrained])
  data.frame(
    Cohort = cohort,
    Group = group,
    N_Proband = model_data$genetic_data$N[1],
    Obs_ConstPTV = observed,
    OE_ConstPTV = observed / expected
  )
}
constrained_rows <- list()
for (cohort in names(cohort_model_names)) {
  for (group_index in seq_along(groups)) {
    model_index <- match(
      cohort_model_names[[cohort]][group_index],
      autism_data$loop_vars$names
    )
    constrained_rows[[length(constrained_rows) + 1L]] <-
      constrained_ptv_row(
        cohort, groups[group_index],
        get_genetic_data(model_index, autism_data)
      )
  }
}

# DDD is a component of the combined autism data but is not fitted as a
# standalone loop entry. Recreate the corresponding count-only specification
# through the same get_genetic_data() path used for Figure 2A.
ddd_data <- autism_data
ddd_index <- length(ddd_data$loop_vars$names) + 1L
combined_index <- match(
  "Combined Probands PTV", ddd_data$loop_vars$names
)
ddd_data$loop_vars$names <- c(
  ddd_data$loop_vars$names, "Combined DDD PTV"
)
ddd_data$loop_vars$datasets[[ddd_index]] <-
  colnames(ddd_data$counts) == "DDD ASD"
ddd_data$loop_vars$N_subset[[ddd_index]] <-
  SummarizedExperiment::colData(ddd_data$counts)["DDD ASD", "N_Proband"]
ddd_data$loop_vars$subsets[[ddd_index]] <-
  ddd_data$loop_vars$subsets[[combined_index]]
ddd_data$loop_vars$mutation_rate[[ddd_index]] <-
  ddd_data$loop_vars$mutation_rate[[combined_index]]
ddd_data$loop_vars$prevalences[[ddd_index]] <-
  ddd_data$loop_vars$prevalences[[combined_index]]
ddd_data$loop_vars$fit_model[[ddd_index]] <- FALSE
constrained_rows[[length(constrained_rows) + 1L]] <-
  constrained_ptv_row(
    "DDD", "All", get_genetic_data(ddd_index, ddd_data)
  )
table2 <- do.call(rbind, constrained_rows)
write_tsv(
  table2,
  "SupplementaryTable2_ConstrainedPTVBurdenByStratum.tsv"
)

# Tables 3 and 8 retain the original model-summary roles, with a corrected,
# symmetric five-column enrichment block for every annotation.
autism_keep <- autism_data$loop_vars$names[
  !grepl("ASC|SPARK|GeneDx|DDID|Non-DDID", autism_data$loop_vars$names)
]
table3 <- model_summary_table(autism_data, BurdenMLE_DN_models_autism, autism_keep)
table8 <- model_summary_table(
  kaplanis_data, BurdenMLE_DN_models_DDD, kaplanis_data$loop_vars$names)
write_tsv(table3, "SupplementaryTable3_AutismModelEstimates.tsv")
write_tsv(
  table8,
  "SupplementaryTable8_DevelopmentalDisorderModelEstimates.tsv"
)

# Table 4: gene-level combined PTV/Mis2 mutational variance.
main_summary <- readRDS(file.path(
  derived_dir, "main_autism", "figure2_summary.rds"))
gnomad <- fread(file.path(
  final_runs_dir, "inputs", "constraint", "gnomad.v2.1.1.lof_metrics.by_gene.txt"),
  select = c("gene", "gene_id"), data.table = FALSE)
table4 <- transform(
  main_summary$gene_mutational_variance,
  gene_name = gnomad$gene[match(gene_id, gnomad$gene_id)]
)[c("gene_id", "gene_name", "mutational_variance",
    "fraction_mutational_variance", "cumulative_fraction")]
names(table4) <- c("gene_ID", "gene_name", "gene_MutVar",
                   "gene_FractionMutVar", "cumulativeFractionMutVar")
write_tsv(
  table4,
  "SupplementaryTable4_GeneLevelMutationalVariance.tsv"
)

fraction_table <- function(x) {
  answer <- x[c(
    "threshold", "ptv", "ptv_lower", "ptv_upper",
    "mis2", "mis2_lower", "mis2_upper"
  )]
  names(answer) <- c(
    "Rate Ratio Threshold", "Fraction of Cases with PTV",
    "Fraction of Cases with PTV Lower95CI",
    "Fraction of Cases with PTV Upper95CI",
    "Fraction of Cases with Mis2",
    "Fraction of Cases with Mis2 Lower95CI",
    "Fraction of Cases with Mis2 Upper95CI"
  )
  answer
}
table5 <- fraction_table(main_summary$fraction_cases)
write_tsv(
  table5,
  "SupplementaryTable5_AutismFractionCasesByRateRatio.tsv"
)

# Table 6: cohort/sex-stratified autism estimates.
cohort_sex <- readRDS(file.path(
  derived_dir, "cohort_ddid", "cohort_sex_summary.rds"))
table6 <- cohort_sex[c(
  "study", "sex", "prevalence", "gamma_scaling_factor",
  "mutvar_ptv", "mutvar_mis2", "mutvar_combined", "mutvar_combined_lower",
  "mutvar_combined_upper", "fraccase_RR5", "fraccase_RR5_lower",
  "fraccase_RR5_upper"
)]
names(table6) <- c(
  "Study", "Sex", "Prevalence", "Scaling_Factor", "MutVar_PTV",
  "MutVar_Mis2", "MutVar_PTV_Mis2", "MutVar_PTV_Mis2_Lower95CI",
  "MutVar_PTV_Mis2_Upper95CI", "Fraction of Cases with RR>5 PTV/Mis2",
  "Fraction of Cases with RR>5 PTV/Mis2 Lower95CI",
  "Fraction of Cases with RR>5 PTV/Mis2 Upper95CI"
)
write_tsv(table6, "SupplementaryTable6_CohortSexEstimates.tsv")

# Table 7: forecasting results. The zero-new-case rows document the shared
# baseline explicitly and are new relative to the original table.
forecast <- fread(file.path(
  derived_dir, "forecasting", "cohort_forecast_summary.tsv"),
  data.table = FALSE)
table7 <- forecast[forecast$threshold == 0, setdiff(names(forecast), "threshold")]
names(table7) <- c(
  "New Data Dataset", "N New Cases", "Total Sample Size",
  "Number FDR Significant Genes (mean)",
  "Number FDR Significant Genes (sd)",
  "Number FDR Significant Genes (lower95CI)",
  "Number FDR Significant Genes (upper95CI)",
  "Number Bonferroni Significant Genes (mean)",
  "Number Bonferroni Significant Genes (sd)",
  "Number Bonferroni Significant Genes (lower95CI)",
  "Number Bonferroni Significant Genes (upper95CI)"
)
write_tsv(table7, "SupplementaryTable7_GeneDiscoveryForecasts.tsv")

# Table 9: DD fraction of cases by rate-ratio threshold.
dd_summary <- readRDS(file.path(
  derived_dir, "autism_dd", "figure4_summary.rds"))
if (is.null(dd_summary$ddd_fraction_by_class)) {
  stop("Autism/DD summary predates class-specific CIs; rerun summarize_autism_dd.R.")
}
table9 <- fraction_table(dd_summary$ddd_fraction_by_class)
write_tsv(
  table9,
  "SupplementaryTable9_DevelopmentalDisorderFractionCasesByRateRatio.tsv"
)

inventory <- data.frame(
  file = c(
    "SupplementaryTable1_SimulationResults.tsv",
    "SupplementaryTable2_ConstrainedPTVBurdenByStratum.tsv",
    "SupplementaryTable3_AutismModelEstimates.tsv",
    "SupplementaryTable4_GeneLevelMutationalVariance.tsv",
    "SupplementaryTable5_AutismFractionCasesByRateRatio.tsv",
    "SupplementaryTable6_CohortSexEstimates.tsv",
    "SupplementaryTable7_GeneDiscoveryForecasts.tsv",
    "SupplementaryTable8_DevelopmentalDisorderModelEstimates.tsv",
    "SupplementaryTable9_DevelopmentalDisorderFractionCasesByRateRatio.tsv"
  ),
  content = c(
    "Simulation estimates (MixSQP)",
    "Constrained-PTV counts and O/E by autism stratum",
    "Autism model estimates",
    "Autism gene-level mutational variance",
    "Autism fraction of cases by RR",
    "Autism cohort/sex estimates",
    "Autism gene-discovery forecasts",
    "Developmental-disorder model estimates",
    "DD fraction of cases by RR"
  ),
  manuscript_status = c(
    "Numbered in original submission",
    "Added during revision",
    rep("Renumbered during revision", 7)
  ),
  manuscript_location = c(
    "Simulation results",
    paste0(
      "First cited after the Figure 2A dataset description; cited again in ",
      "the cohort/sex heterogeneity paragraph"
    ),
    "Autism mutational-variance and effective-penetrance results",
    "Autism polygenicity results",
    "Autism fraction-of-cases results",
    "Autism cohort/sex heterogeneity results",
    "Future autism gene-discovery results",
    "Developmental-disorder comparison",
    "Developmental-disorder comparison"
  )
)
write_tsv(inventory, "supplementary_table_inventory.tsv")

cat("Supplementary tables written to", output_dir, "\n")
print(data.frame(
  table = basename(inventory$file),
  rows = c(nrow(table1), nrow(table2), nrow(table3), nrow(table4),
           nrow(table5), nrow(table6), nrow(table7), nrow(table8),
           nrow(table9))
), row.names = FALSE)
