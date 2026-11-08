# Fresh simulation batches; original scientific functions, independent batch RNG.
root <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT"), mustWork = TRUE)
repo <- normalizePath(Sys.getenv("BURDENMLEDN_REPO_ROOT"))
args <- commandArgs(TRUE)
stage <- args[1]
run_date <- Sys.getenv("RUN_DATE")
seed <- as.integer(Sys.getenv("SIMULATION_SEED", "20260715"))
out <- file.path(root, "outputs/data")
jobs <- file.path(root, "jobs/simulations")
dir.create(jobs, recursive = TRUE, showWarnings = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
library(pbapply)
for (f in list.files(file.path(repo, "R"), "[.]R$", full.names = TRUE)) source(f)
source(file.path(repo, "analysis/scripts/secondary_analysis_functions.R"))
source(file.path(repo, "analysis/scripts/model_jobs.R"))
source(file.path(repo, "analysis/scripts/simulation_functions.R"))
distributions <- c("null", "half_uniform", "half_infinitesimal", "point_normal", "oligogenic")
sample_sizes <- c(1000, 5000, 15000, 38680, 50000)
inventory <- expand.grid(batch = 1:20, N = sample_sizes, distribution = distributions,
                         stringsAsFactors = FALSE)
inventory$id <- seq_len(nrow(inventory))
inventory$seed <- seed + inventory$id
inventory$first <- (inventory$batch - 1L) * 25L + 1L
inventory$last <- inventory$batch * 25L
inventory$compare_EM <- inventory$last <= 100L
cache <- file.path(jobs, "inputs.rds")
if (stage == "prepare") {
  # Same empirical reference as get_genetic_data(5, autism_data) in the serial study.
  x <- readRDS(file.path(root, file.path("jobs", run_date, "main", "Combined_Probands_PTV_cap100_p01_m005.rds")))
  saveRDS(list(data = x$arguments$input_data, features = x$arguments$features,
               prevalence = x$arguments$prevalence), cache)
  write.table(inventory, file.path(jobs, "inventory.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
  set.seed(seed)
  params <- list(null=list(),
    half_uniform=list(endpoints=seq(0,log(100),length.out=10),weights=draw_dirichlet(rep(1,10))),
    half_infinitesimal=list(sd=runif(1,.5,1)),
    point_normal=list(pi0=runif(1,.85,.95),mean=runif(1,1.5,2.5),sd=runif(1,.2,.4)),
    oligogenic=list(prop_causal=runif(1,.001,.01),effect_size=runif(1,log(5),log(50))))
  effect_size_df <- do.call(rbind, lapply(distributions, function(d)
    data.frame(effect_size=generate_effect_sizes(nrow(x$arguments$input_data),d,params[[d]]), distribution=d)))
  save(effect_size_df, file=file.path(out,sprintf("effect_size_samples_%s_seed%s.RData",run_date,seed)))
} else if (stage == "worker") {
  id <- as.integer(args[2]); stopifnot(id %in% inventory$id)
  job <- inventory[id, ]; x <- readRDS(cache)
  params <- switch(job$distribution,
    null=list(),
    half_uniform=list(endpoints=seq(0,log(100),length.out=10),weights=function() draw_dirichlet(rep(1,10))),
    half_infinitesimal=list(sd=function() runif(1,.5,1)),
    point_normal=list(pi0=function() runif(1,.85,.95),mean=function() runif(1,1.5,2.5),sd=function() runif(1,.2,.4)),
    oligogenic=list(prop_causal=function() runif(1,.001,.01),effect_size=function() runif(1,log(5),log(50))))
  params <- c(params,list(prevalence=x$prevalence,RR_threshold=5))
  path <- file.path(jobs,sprintf("batch_%03d.rds",id))
  if (file.exists(path)) stop("Refusing to overwrite completed batch: ",path)
  z <- list()
  for (method in c("mixsqp", if(job$compare_EM) "EM")) {
    set.seed(job$seed)
    z[[method]] <- run_simulation(job$distribution,params,x$data,x$features,
                                  job$N,n_sims=25,optimizer=method)
    for(i in seq_along(z[[method]])) z[[method]][[i]]$simulation_index_within_N <- job$first+i-1L
  }
  if(job$compare_EM) for(i in 1:25) {
    stopifnot(identical(z$mixsqp[[i]]$true_effects,z$EM[[i]]$true_effects),
              identical(z$mixsqp[[i]]$simulation_data_hash,z$EM[[i]]$simulation_data_hash),
              identical(z$mixsqp[[i]]$true_frac_cases,z$EM[[i]]$true_frac_cases))
  }
  atomic_job_rds(list(job=job,results=z),path)
} else if (stage == "collect") {
  all <- lapply(inventory$id,function(id) {
    z <- readRDS(file.path(jobs,sprintf("batch_%03d.rds",id)))
    stopifnot(identical(z$job,inventory[id, ])); z
  })
  for(method in c("mixsqp","EM")) {
    results <- setNames(lapply(distributions,function(d) {
      blocks <- Filter(function(z) z$job$distribution==d && !is.null(z$results[[method]]),all)
      unlist(lapply(blocks,function(z) z$results[[method]]),recursive=FALSE)
    }),distributions)
    expected <- if(method=="mixsqp") 500L else 100L
    for(d in distributions) for(n in sample_sizes) {
      z <- Filter(function(x) x$N==n,results[[d]])
      stopifnot(length(z)==expected,
                identical(vapply(z,`[[`,integer(1),"simulation_index_within_N"),seq_len(expected)),
                all(vapply(z,function(x) isTRUE(x$optimizer_converged),logical(1))))
    }
    simulation_metadata <- list(optimizer=method,rng_seed=seed,run_date=run_date,
      sample_sizes=sample_sizes,distributions=distributions,
      simulations_per_distribution_and_N=expected,optimizer_rng_isolated=TRUE,
      batch_seeds=inventory$seed,batch_size=25L,comparison_replicates=100L,
      component_endpoints=seq(0,log(100),length.out=10))
    save(results,simulation_metadata,file=file.path(out,sprintf("simulation_results_%s_%s_seed%s.RData",tolower(method),run_date,seed)))
  }
} else stop("Expected prepare, worker, or collect")
