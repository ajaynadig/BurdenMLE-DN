root <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT"))
repo <- normalizePath(Sys.getenv("BURDENMLEDN_REPO_ROOT"))
dir.create(file.path(root,"provenance"),showWarnings=FALSE)
if (!dir.exists(file.path(root,"inputs"))) stop("Provide authorized inputs using analysis/input_manifest.tsv")
packages <- c("pbapply", "data.table", "ggplot2", "patchwork", "SummarizedExperiment",
              "readr", "dplyr", "tidyr", "scales", "openxlsx", "digest", "MetBrewer")
missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop("Install analysis dependencies before running: ",paste(missing,collapse=", "))
lib <- file.path(root,"R-library")
status <- system2(file.path(R.home("bin"),"R"),c("CMD","INSTALL","--library",shQuote(lib),shQuote(repo)))
if(status != 0) stop("Package installation failed")
writeLines(c(capture.output(sessionInfo()),capture.output(.libPaths())),file.path(root,"provenance/R_environment.txt"))
write.table(as.data.frame(installed.packages()[,c("Package","Version","LibPath")]),
            file.path(root,"provenance/R_packages.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
source(file.path(repo,"analysis/workflow/record_inputs.R"))
draws <- Sys.getenv("BURDENMLEDN_PREVIOUS_BOOTSTRAPS")
if(nzchar(draws)) {
  paths <- file.path(draws,paste0(c("main_autism","main_ddd","no_ces","no_overlap","age"),".rds"))
  stopifnot(all(file.exists(paths)))
  write.table(data.frame(path=normalizePath(paths),sha256=vapply(paths,digest::digest,character(1),algo="sha256",file=TRUE)),
              file.path(root,"provenance/bootstrap_inputs.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
}
# The CDC archive uses Deflate64. Extract once with a compatible unzip before
# launching the study; a compute node need not install another decompressor.
archive <- file.path(root,"inputs/demographics/Nat2024us.zip")
if(!file.exists(archive)) stop("Provide the public CDC natality archive and extracted text (see analysis README).")
members <- unzip(archive,list=TRUE)
i <- which.max(members$Length)
text <- file.path(dirname(archive),members$Name[i])
if(!file.exists(text) || file.info(text)$size != members$Length[i]) stop("Extract the complete natality archive before launch.")
