args <- commandArgs(trailingOnly = TRUE)
library(pbapply)
if (length(args) != 1L) stop("Usage: Rscript run_model_job.R /path/to/job.rds")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo_dir <- dirname(dirname(dirname(normalizePath(script))))
for (f in c("BurdenMLE_DN.R", "EM.R", "MixSQP.R", "estimate_mutvar.R",
            "io.R", "likelihoods.R", "model.R")) source(file.path(repo_dir, "R", f))
source(file.path(repo_dir, "analysis/scripts/model_jobs.R"))
execute_model_job(args[1], repo_dir)
