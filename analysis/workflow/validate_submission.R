root <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT"))
cohort <- readRDS(file.path(root,"outputs/derived/cohort_ddid/cohort_sex_summary.rds"))
stopifnot(all(cohort$prevalence[cohort$sex=="Both" & cohort$study %in% c("ASC","SPARK")]==.0276),
          all(cohort$prevalence[cohort$sex=="Both" & cohort$study %in% c("GeneDx","GeneDX")]==.01))
strata <- readRDS(file.path(root,"outputs/derived/cohort_ddid/ddid_summary.rds"))
# These summaries are generated from the selected fits; assert_fit_prevalence
# checks their prevalence at generation rather than relying on a prior run.
dest <- file.path(root,"submission")
d <- read.delim(file.path(dest,"ManuscriptEstimateChecklist.tsv"))
required <- c(R2_cohort_summary=9L,R2_GeneDx_reporting_correction=8L,R2_paternal_age_summary=1L,
              R2_500_replicate_summary=20L,R2_CES_excluded_sensitivity=80L)
for(s in names(required)) stopifnot(sum(d$update_status==s)==required[[s]])
stopifnot(nrow(d)==214L,!anyDuplicated(d[c("manuscript_section","paragraph_anchor","estimand","stratum","source_key")]),
          setequal(d$estimate[d$estimand=="Maximum fitted rate ratio"],c(100,250)))
files <- list.files(dest,recursive=TRUE,full.names=TRUE)
files <- files[!file.info(files)$isdir]
stopifnot(all(file.info(files)$size>0))
write.table(data.frame(file=substring(files,nchar(root)+2),sha256=vapply(files,digest::digest,character(1),algo="sha256",file=TRUE)),
            file.path(root,"provenance/output_manifest.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
cat("Submission completeness and prevalence checks passed. Visual review remains required before release.\n")
