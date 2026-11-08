if (dir.exists("analysis/workflow")) {
root <- normalizePath(".")
repo <- root
for (f in c("likelihoods.R","io.R","model.R")) source(file.path(repo,"R",f))
# Exercise the actual forecasting functions, not a simplified RNG surrogate.
for (e in parse(file.path(repo,"analysis/scripts/forecasting_script_revision.R"))) {
  if (is.call(e) && identical(e[[1]],as.name("<-")) && is.call(e[[3]]) &&
      identical(e[[3]][[1]],as.name("function"))) eval(e)
}
data <- read.csv(file.path(repo,"inst/extdata/synthetic_example.csv"),row.names="gene")
data <- data[seq_len(60),c("case_count","case_rate","N")]
data$expected_count <- 2 * data$N * data$case_rate
fit <- initialize_model(poisson_uniform_likelihood, data, c(0,log(20)), features=NULL,grid_size=10L)
fit$delta[,] <- c(.7,.3)
fit$posterior_gene_estimates <- data.frame(Case_Count=data$case_count,
  Expected_Count=data$expected_count,row.names=rownames(data))
class(fit) <- c("BurdenMLEDN_fit","list")
sampler <- posterior_gene_sampler(fit,data)
scenario_list <- list(A=list(data=data,sampler=sampler,scaling=1,group="cohort"),
                      B=list(data=data,sampler=sampler,scaling=.8,group="ddid"))
tasks <- expand.grid(n_new_case=c(0,1000,5000),iteration=1:4,
                     scenario=names(scenario_list),stringsAsFactors=FALSE)
seed <- 24312342L
serial <- lapply(seq_len(nrow(tasks)),run_task)
for (cores in c(2L,4L)) {
  parallel <- parallel::mclapply(seq_len(nrow(tasks)),run_task,
                                mc.cores=cores,mc.set.seed=FALSE)
  reverse <- parallel::mclapply(rev(seq_len(nrow(tasks))),run_task,
                               mc.cores=cores,mc.set.seed=FALSE)
  stopifnot(identical(serial,parallel),identical(serial,rev(reverse)),
    identical(summarize_forecast(do.call(rbind,serial)),
              summarize_forecast(do.call(rbind,parallel))))
}
retry <- lapply(c(7L,1L,24L),run_task)
stopifnot(identical(retry,serial[c(7L,1L,24L)]))
cat("Forecast task outputs and summaries identical: 1, 2, 4 cores; reverse scheduling; retries.\n")

}
