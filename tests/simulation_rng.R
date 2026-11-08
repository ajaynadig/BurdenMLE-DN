if (dir.exists("analysis/workflow")) {
root <- normalizePath(".")
repo <- root
library(pbapply)
for(f in list.files(file.path(repo,"R"),"[.]R$",full.names=TRUE)) source(f)
source(file.path(repo,"analysis/scripts/secondary_analysis_functions.R"))
source(file.path(repo,"analysis/scripts/simulation_functions.R"))
dat <- read.csv(file.path(repo,"inst/extdata/synthetic_example.csv"),row.names="gene")
features <- matrix(1,nrow(dat),1,dimnames=list(rownames(dat),"all"))
params <- list(pi0=.9,mean=2,sd=.3,prevalence=.022,RR_threshold=5)
set.seed(20260716)
a <- run_simulation("point_normal",params,dat,features,5000,n_sims=4,optimizer="mixsqp")
set.seed(20260716)
b <- c(run_simulation("point_normal",params,dat,features,5000,n_sims=2,optimizer="mixsqp"),
       run_simulation("point_normal",params,dat,features,5000,n_sims=2,optimizer="mixsqp"))
set.seed(20260716)
c <- run_simulation("point_normal",params,dat,features,5000,n_sims=4,optimizer="EM")
for(i in 1:4) {
  stopifnot(identical(a[[i]]$simulation_data_hash,b[[i]]$simulation_data_hash),
            identical(a[[i]]$simulation_data_hash,c[[i]]$simulation_data_hash),
            identical(a[[i]]$true_effects,c[[i]]$true_effects),
            identical(a[[i]]$estimated_mutvar,b[[i]]$estimated_mutvar))
}
cat("Simulation split/replay and EM pairing passed.\n")

}
