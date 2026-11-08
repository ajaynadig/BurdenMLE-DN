# Rebuild submission exports from this run only; no R1/R2 result reuse.
root <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT"))
tables <- file.path(root,"outputs/tables/supplementary")
dest <- file.path(root,"submission")
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
table_dest <- file.path(dest,"Supplementary Tables")
dir.create(table_dest,showWarnings=FALSE)
files <- list.files(tables,"^SupplementaryTable[1-9]_.*[.]tsv$",full.names=TRUE)
stopifnot(length(files)==9)
old <- as.integer(sub("SupplementaryTable([0-9]+)_.*","\\1",basename(files)))
new <- ifelse(old==1,1,old+1)
# Map each table explicitly; the bias/precision summary takes position 2.
names_new <- vapply(seq_along(files),function(i)
  sub("^SupplementaryTable[0-9]+_",paste0("SupplementaryTable",new[i],"_"),basename(files[i])),character(1))
stopifnot(all(file.copy(files,file.path(table_dest,names_new),overwrite=TRUE)))
x <- read.delim(file.path(tables,"SupplementaryTable1_SimulationResults.tsv"))
stopifnot(nrow(x)==10000,all(is.finite(x$True.Mutational.Variance)),all(x$True.Mutational.Variance>0))
groups <- split(x,interaction(x$Distribution,x$N,drop=TRUE))
precision <- do.call(rbind,lapply(groups,function(z) {
  stopifnot(nrow(z)==500L)
  ratio <- z$Estimated.Mutational.Variances/z$True.Mutational.Variance
  stopifnot(all(is.finite(ratio)))
  q <- quantile(ratio,c(.05,.95))
  data.frame(Architecture=z$Distribution[1],N=z$N[1],Replicates=nrow(z),
    Mean_Estimate_Truth_Ratio=mean(ratio),Mean_Relative_Error=mean(ratio-1),
    Relative_RMSE=sqrt(mean((ratio-1)^2)),Percentile_5=q[1],Percentile_95=q[2],
    Percentile_95_minus_5=diff(q))
}))
precision_name <- "SupplementaryTable2_SimulationBiasPrecision.tsv"
write.table(precision,file.path(table_dest,precision_name),sep="\t",row.names=FALSE,quote=FALSE)
map <- rbind(data.frame(old=old,new=new,file=names_new),
             data.frame(old=NA_integer_,new=2,file=precision_name))
map <- map[order(map$new),]
write.table(map,file.path(dest,"TableNumbering.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
# Read back the workbook and compare every table to its TSV source.
wb <- openxlsx::createWorkbook()
for(i in seq_len(nrow(map))) {
  d <- read.delim(file.path(table_dest,map$file[i]),check.names=FALSE)
  sheet <- paste("Table",map$new[i]); openxlsx::addWorksheet(wb,sheet)
  openxlsx::writeData(wb,sheet,d,keepNA=TRUE)
  openxlsx::freezePane(wb,sheet,firstRow=TRUE)
  openxlsx::setColWidths(wb,sheet,cols=seq_len(ncol(d)),widths=18)
}
workbook <- file.path(table_dest,"SupplementaryTables.xlsx")
openxlsx::saveWorkbook(wb,workbook,overwrite=TRUE)
for(i in seq_len(nrow(map))) {
  expected <- read.delim(file.path(table_dest,map$file[i]),check.names=FALSE)
  actual <- openxlsx::read.xlsx(workbook,sheet=i,check.names=FALSE,skipEmptyRows=FALSE,skipEmptyCols=FALSE)
  stopifnot(isTRUE(all.equal(expected,actual,check.attributes=FALSE,tolerance=1e-10)))
}
for(kind in c("main","supplementary")) {
  target <- file.path(dest,paste0(kind,"_figures")); dir.create(target,showWarnings=FALSE)
  source <- list.files(file.path(root,"outputs/figures",kind),"[.]pdf$",full.names=TRUE)
  stopifnot(length(source)>0,all(file.copy(source,target,overwrite=TRUE)))
}
stopifnot(file.copy(file.path(root,"outputs/derived/upper_rr_no_CES/UpperRRSweep_all_metrics.pdf"),
  file.path(dest,"supplementary_figures/SupplementaryFigure17_UpperRRBounds.pdf"),overwrite=TRUE))
inventory <- list.files(file.path(dest,"supplementary_figures"),"^SupplementaryFigure[0-9]+_.*[.]pdf$")
numbers <- as.integer(sub("SupplementaryFigure([0-9]+)_.*","\\1",inventory))
stopifnot(identical(sort(numbers),1:17))
write.table(data.frame(figure=numbers,file=inventory),file.path(dest,"FigureNumbering.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
estimates <- list.files(file.path(root,"outputs"),"^ManuscriptEstimateChecklist[.]tsv$",recursive=TRUE,full.names=TRUE)
stopifnot(length(estimates)==1L)
d <- read.delim(estimates,check.names=FALSE)
for(col in names(d)) if(is.character(d[[col]])) {
  for(i in 9:2) d[[col]] <- gsub(paste0("Supplementary Table ",i,"\\b"),paste0("Supplementary Table ",i+1),d[[col]],perl=TRUE)
}
write.table(d,file.path(dest,"ManuscriptEstimateChecklist.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
writeLines(c("# Fresh R2 manuscript handoff", "",
  "Generated from the fresh run. Figures, tables and workbook still require visual inspection before freezing.",
  "", "1. Use ManuscriptEstimateChecklist.tsv for numerical manuscript edits.",
  "2. TableNumbering.tsv maps old Tables 2–9 to 3–10; the simulation precision summary is Table 2.",
  "3. Supplementary Figure 6 is the two-panel TADA density plot. Figure 17 is the CES-excluded upper-bound sensitivity.",
  "4. Simulation precision statistics now use 500 replicates per scenario; replace all prior quoted simulation summaries.",
  "5. ASC/SPARK use population prevalence; GeneDx uses clinical prevalence. Check the fresh DDID/cohort summaries and paternal-age comparison.",
  "6. The public Kaplanis example is a separate demonstration, not the manuscript DD data."),file.path(dest,"READ_ME_FIRST.md"))
