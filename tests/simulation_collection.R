if (dir.exists("analysis/workflow")) {
# Synthetic file-format fixture only; never creates scientific run results.
main <- function() {
  root <- normalizePath(".")
  tmp <- tempfile("simulation-collection-"); dir.create(tmp)
  on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
  dir.create(file.path(tmp,"code"))
  stopifnot(file.symlink(root,file.path(tmp,"code/BurdenMLE-DN")))
  dir.create(file.path(tmp,"jobs/simulations"),recursive=TRUE)
  distributions <- c("null","half_uniform","half_infinitesimal","point_normal","oligogenic")
  inventory <- expand.grid(batch=1:20,N=c(1000,5000,15000,38680,50000),distribution=distributions,stringsAsFactors=FALSE)
  inventory$id <- seq_len(nrow(inventory)); inventory$seed <- 20260715L+inventory$id
  inventory$first <- (inventory$batch-1L)*25L+1L; inventory$last <- inventory$batch*25L
  inventory$compare_EM <- inventory$last<=100L
  for(id in inventory$id) {
    job <- inventory[id,]; results <- list()
    for(method in c("mixsqp",if(job$compare_EM)"EM")) results[[method]] <- lapply(job$first:job$last,function(i)
      list(N=job$N,simulation_index_within_N=as.integer(i),optimizer_converged=TRUE,
           true_mutvar=i/10000,estimated_mutvar=i/10000,true_frac_cases=.1,estimated_frac_cases=.1,
           true_peneff=.2,estimated_peneff=.2,optimizer=method,active_components_by_stratum=1))
    saveRDS(list(job=job,results=results),file.path(tmp,"jobs/simulations",sprintf("batch_%03d.rds",id)))
  }
  old <- Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT")
  Sys.setenv(BURDENMLEDN_ANALYSIS_ROOT=tmp, BURDENMLEDN_REPO_ROOT=root, RUN_DATE="fixture")
  status <- system2(file.path(R.home("bin"),"Rscript"),c("--vanilla",shQuote(file.path(root,"analysis/workflow/simulation_batches.R")),"collect"))
  Sys.setenv(BURDENMLEDN_ANALYSIS_ROOT=old)
  stopifnot(status==0)
  env <- new.env(parent=globalenv())
  for(e in parse(file.path(root,"analysis/scripts/analyze_simulations.R")))
    if(is.call(e)&&identical(e[[1]],as.name("<-"))&&identical(e[[2]],as.name("create_summary_df"))) eval(e,env)
  a <- new.env(); b <- new.env()
  load(file.path(tmp,"outputs/data/simulation_results_mixsqp_fixture_seed20260715.RData"),a)
  load(file.path(tmp,"outputs/data/simulation_results_em_fixture_seed20260715.RData"),b)
  env$summary_df <- env$create_summary_df(a$results,"mixsqp")
  env$comparison_summary <- env$create_summary_df(b$results,"EM")
  env$optimizer_label <- "mixsqp"; env$compare_optimizer <- "EM"
  env$get_arg <- function(...) "100"
  lines <- readLines(file.path(root,"analysis/scripts/analyze_simulations.R"))
  start <- grep('^  pair_columns <-',lines); end <- grep('^  primary_suffix <-',lines)-1L
  eval(parse(text=lines[start:end]),env)
  stopifnot(nrow(env$summary_df)==10000L,nrow(env$paired_summary)==2000L)
  cat("500-replicate collection and exactly 100 matched comparator replicates per non-null scenario passed.\n")
}
main()

}
