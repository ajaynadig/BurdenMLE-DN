if (dir.exists("analysis/workflow")) {
# Test export mechanics with explicitly synthetic fixtures, never scientific outputs.
main <- function() {
  root <- normalizePath(".")
  tmp <- tempfile("submission-export-"); dir.create(tmp)
  on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
  tables <- file.path(tmp,"outputs/tables/supplementary"); dir.create(tables,recursive=TRUE)
  x <- expand.grid(Distribution=c("half_uniform","half_infinitesimal","point_normal","oligogenic"),N=c(1000,5000,15000,38680,50000),replicate=1:500)
  x$True.Mutational.Variance <- .01; x$Estimated.Mutational.Variances <- .011
  write.table(x,file.path(tables,"SupplementaryTable1_SimulationResults.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  for(i in 2:9) write.table(data.frame(Gene=c("A","B"),Estimate=c(.01,.02),CI=c(NA,.03)),
     file.path(tables,paste0("SupplementaryTable",i,"_Fixture.tsv")),sep="\t",row.names=FALSE,quote=FALSE)
  for(kind in c("main","supplementary")) {
    path <- file.path(tmp,"outputs/figures",kind); dir.create(path,recursive=TRUE)
    for(i in if(kind=="main") 1:4 else 1:16)
      writeLines("TEST FIXTURE NOT A FIGURE",file.path(path,paste0(if(kind=="main")"Figure" else "SupplementaryFigure",i,"_Fixture.pdf")))
  }
  path <- file.path(tmp,"outputs/derived/upper_rr_no_CES"); dir.create(path,recursive=TRUE)
  writeLines("TEST FIXTURE NOT A FIGURE",file.path(path,"UpperRRSweep_all_metrics.pdf"))
  path <- file.path(tmp,"outputs/tables/manuscript_values"); dir.create(path,recursive=TRUE)
  write.table(data.frame(text="Supplementary Table 2; Supplementary Table 9"),
    file.path(path,"ManuscriptEstimateChecklist.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  old <- Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT"); Sys.setenv(BURDENMLEDN_ANALYSIS_ROOT=tmp, BURDENMLEDN_REPO_ROOT=root, RUN_DATE="fixture")
  status <- system2(file.path(R.home("bin"),"Rscript"),c("--vanilla",shQuote(file.path(root,"analysis/workflow/assemble_submission.R"))))
  Sys.setenv(BURDENMLEDN_ANALYSIS_ROOT=old)
  stopifnot(status==0)
  z <- read.delim(file.path(tmp,"submission/ManuscriptEstimateChecklist.tsv"))
  stopifnot(z$text=="Supplementary Table 3; Supplementary Table 10")
  cat("Submission numbering, 500-replicate summary, and workbook round-trip fixture passed.\n")
}
main()

}
