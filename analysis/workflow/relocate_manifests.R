root <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT"))
repo <- normalizePath(Sys.getenv("BURDENMLEDN_REPO_ROOT"))
source(file.path(repo,"analysis/scripts/model_artifacts.R"))
for(mode in c("main","no_ces","no_overlap","age")) {
  path <- file.path(root,"outputs/data",paste0("model_manifest_",mode,"_",Sys.getenv("RUN_DATE"),".rds"))
  m <- read_model_manifest(path)
  for(name in names(m$artifacts)) {
    m$artifacts[[name]]$path <- normalizePath(file.path(root,"outputs/data",basename(m$artifacts[[name]]$path)),mustWork=TRUE)
  }
  saveRDS(m,path)
  for(name in names(m$artifacts)) load_model_artifact(path,name,envir=new.env())
}
