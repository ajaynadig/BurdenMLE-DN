# Sourced by make_manuscript_estimate_table.R after the primary checklist.
# All values come from this run's derived products; no historical values are copied.
checklist$stratum[checklist$estimand == "Maximum fitted rate ratio"] <- "Primary autism analysis"
add_row("Methods", "Implementation details", "Maximum fitted rate ratio",
        "Primary developmental-disorder analysis", 250, unit = "rate ratio",
        manuscript_ready = "250", source_file = file.path(repo_dir, "analysis/scripts/set_up_kaplanis.R"),
        source_key = "DD model specification")
cohort_file <- file.path(cohort_ddid_dir, "cohort_sex_summary.tsv")
cohort <- read_tsv(cohort_file)
stopifnot(nrow(cohort) == 9L)
for (i in seq_len(nrow(cohort))) {
  z <- cohort[i, ]
  add_row("Results", "Cohort/sex heterogeneity", "PTV + Mis2 mutational variance",
          paste(z$study, z$sex), z$mutvar_combined,
          z$mutvar_combined_lower, z$mutvar_combined_upper, unit = "proportion",
          manuscript_ready = format_percent_ci(z$mutvar_combined, z$mutvar_combined_lower, z$mutvar_combined_upper),
          status = "R2_cohort_summary", source_file = cohort_file, source_key = paste(z$study, z$sex),
          notes = paste("Prevalence", z$prevalence, "; fresh fit with historical bootstrap draws."))
}
ddid_file <- file.path(cohort_ddid_dir, "ddid_summary.tsv")
ddid <- subset(read_tsv(ddid_file), Dataset == "GeneDx")
stopifnot(nrow(ddid) == 4L)
for (i in seq_len(nrow(ddid))) {
  z <- ddid[i, ]
  for (metric in c("PTV_peneff", "PTV_effRR")) {
    pen <- metric == "PTV_peneff"
    add_row("Results", "Maximally stratified DD/ID analyses",
            if (pen) "PTV mutation-weighted effective penetrance" else "PTV effective rate ratio",
            z$stratum, z[[metric]], z[[paste0(metric,"_lower")]], z[[paste0(metric,"_upper")]],
            unit = if (pen) "proportion" else "rate ratio",
            manuscript_ready = if (pen) format_percent_ci(z[[metric]], z[[paste0(metric,"_lower")]], z[[paste0(metric,"_upper")]]) else
              format_number_ci(z[[metric]], z[[paste0(metric,"_lower")]], z[[paste0(metric,"_upper")]], 2),
            status = "R2_GeneDx_reporting_correction", source_file = ddid_file, source_key = z$stratum,
            notes = "Fresh fit uses the intended GeneDx stratum prevalence; Supplementary Figure 12.")
  }
}
add_row("Results", "Population and clinical genetic characterization",
        "SPARK paternal-age scalar sensitivity", "PTV + Mis2; prevalence 2.76%",
        paternal_combined$paternal_age_adjusted_mutvar, unit = "proportion",
        manuscript_ready = paste0(format_plain(100*paternal_combined$reported_mutvar,3), "% to ",
                                  format_plain(100*paternal_combined$paternal_age_adjusted_mutvar,3), "%"),
        status = "R2_paternal_age_summary", source_file = paternal_bias_file,
        source_key = "Combined PTV + Mis2", notes = "Scalar sensitivity, not a mutation-rate refit; PTV indel scaling included.")

simulation_file <- file.path(tables_dir,"supplementary","SupplementaryTable1_SimulationResults.tsv")
simulation <- read_tsv(simulation_file)
names(simulation) <- make.names(names(simulation))
groups <- split(simulation, interaction(simulation$Distribution, simulation$N, drop=TRUE))
stopifnot(length(groups)==20L)
for (z in groups) {
  ratio <- z$Estimated.Mutational.Variances / z$True.Mutational.Variance
  stopifnot(nrow(z)==500L, all(is.finite(ratio)), all(z$True.Mutational.Variance>0))
  q <- quantile(ratio,c(.05,.95))
  label <- switch(as.character(z$Distribution[1]), point_normal="Point-normal",
                  half_infinitesimal="Half-infinitesimal", oligogenic="Oligogenic", half_uniform="Half-uniform")
  add_row("Results", "Simulation results", "Mean relative mutational-variance error",
          paste0(label,"; N = ", format(z$N[1],big.mark=",",scientific=FALSE,trim=TRUE),"; full simulated range"),
          mean(ratio-1), unit="relative bias",
          manuscript_ready=paste0(format_plain(100*mean(ratio-1),1),"%; empirical ratio range ",
                                 format_plain(q[1],2),"–",format_plain(q[2],2)),
          status="R2_500_replicate_summary", source_file=simulation_file,
          source_key=paste0(z$Distribution[1],"; N=",z$N[1]),
          notes="500 replicates. Signed mean relative error; 5th–95th percentile ratio range is replicate dispersion, not a CI. See the simulation bias/precision table in submission exports.")
}

# The accepted sensitivity is the CES-excluded sweep, replacing the preliminary
# all-gene broad/capped experiment. Keep that distinction explicit in the checklist.
sweep_file <- file.path(derived_dir,"upper_rr_no_CES","full_metric_summary.tsv")
sweep <- read_tsv(sweep_file)
sweep <- subset(sweep, metric != "loglik")
stopifnot(nrow(sweep)==60L, !anyDuplicated(sweep[c("dataset","RR_cap","metric")]))
metric_names <- c(mutvar="Mutational variance", effective_penetrance="Mutation-weighted effective penetrance",
                  fraction_cases_RR5="Fraction of cases with RR > 5")
for (i in seq_len(nrow(sweep))) {
  z <- sweep[i, ]
  add_row("Sensitivity", "CES-excluded upper-bound sensitivity", metric_names[[z$metric]],
          paste(z$dataset,paste0("cap=",format(z$RR_cap,digits=8))),
          z$estimate,z$lower95,z$upper95,unit="proportion",
          manuscript_ready=format_percent_ci(z$estimate,z$lower95,z$upper95),
          status="R2_CES_excluded_sensitivity", source_file=sweep_file,
          source_key=paste(z$dataset,z$RR_cap,z$metric),
          notes="Supplementary Figure 17; CES genes excluded. Variance/case fractions combine PTV+Mis2 for both ASD and DD; penetrance is PTV only. Not a primary estimate.")
}
ll_file <- file.path(derived_dir,"upper_rr_no_CES","pooled_estimates.tsv")
ll <- read_tsv(ll_file)
stopifnot(nrow(ll)==20L, setequal(unique(ll$dataset),c("ASD: PTV + Mis2","DD: PTV + Mis2")))
for (i in seq_len(nrow(ll))) {
  z <- ll[i, ]
  add_row("Sensitivity", "CES-excluded upper-bound sensitivity", "Log likelihood relative to best bound",
          paste(z$dataset,paste0("cap=",format(z$RR_cap,digits=8))),z$loglik_difference_from_best,
          unit="log likelihood",manuscript_ready=format_plain(z$loglik_difference_from_best,2),
          source_file=ll_file,source_key=paste(z$dataset,z$RR_cap),status="R2_CES_excluded_sensitivity")
}
checklist$order <- seq_len(nrow(checklist))
key <- c("manuscript_section","paragraph_anchor","estimand","stratum","source_key")
stopifnot(!anyDuplicated(checklist[key]),
          sum(checklist$update_status=="R2_cohort_summary")==9L,
          sum(checklist$update_status=="R2_GeneDx_reporting_correction")==8L,
          sum(checklist$update_status=="R2_paternal_age_summary")==1L,
          sum(checklist$update_status=="R2_500_replicate_summary")==20L,
          sum(checklist$update_status=="R2_CES_excluded_sensitivity")==80L)
