# Fresh CES-excluded upper-bound fits, one model/cap per worker.
root <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT"))
repo <- normalizePath(Sys.getenv("BURDENMLEDN_REPO_ROOT"))
args <- commandArgs(TRUE); stage <- args[1]
out <- file.path(root,"outputs/derived/upper_rr_no_CES")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
for(f in list.files(file.path(repo,"R"),"[.]R$",full.names=TRUE)) source(f)
source(file.path(repo,"analysis/scripts/secondary_analysis_functions.R"))
source(file.path(repo,"analysis/scripts/model_jobs.R"))
source(file.path(repo,"analysis/workflow/upper_bound_helpers.R"))
moments <- function(m, d, legacy_grid = FALSE) {
  u <- if (legacy_grid) seq(.05, 1, by = 1 / m$grid_size) else midpoint_grid(m$grid_size)
  a <- b <- f <- matrix(0, nrow(d), length(m$component_endpoints))
  for (j in seq_along(m$component_endpoints)) {
    rr <- exp(u * m$component_endpoints[j])
    w <- normalize_log_likelihood_rows(dpois(
      matrix(d$case_count, nrow(d), length(u)), outer(d$expected_count, rr), log = TRUE
    ))
    a[, j] <- drop(w %*% (rr * (rr - 1)))
    b[, j] <- drop(w %*% (rr - 1))
    f[, j] <- drop(w %*% (rr * (rr > 5)))
  }
  list(numerator = a, denominator = b, fraction = f)
}
measure <- function(m, d, q, prevalence, scale) {
  cp <- component_posterior_probabilities(m)
  c(mutvar = estimate_mutvar_trio(m, d, prevalence)$total_mutvar * scale,
    effective_penetrance = effective_penetrance_from_moments(m, q, prevalence, d$case_rate),
    fraction_cases_RR5 = sum(2 * d$case_rate * rowSums(cp * q$fraction)) * scale)
}
subset_moments <- function(q, idx) lapply(q, function(a) a[idx, , drop = FALSE])


workspace <- file.path(out,"workspace_private.rds")
if(stage=="prepare") {
  source(file.path(repo,"analysis/scripts/set_up_asc.R"))
  asd_scale <- autism_data$ptv_scale_factor
  dd <- read.table(file.path(root,"inputs/developmental_disorder/kaplanis_variants_GRCh38_annotated_formatted_2024-11-23.txt"),header=TRUE)
  dd_scale <- sum(dd$isPTV,na.rm=TRUE)/sum(dd$isPTV & !dd$isIndel,na.rm=TRUE)
  paths <- c("Combined_Probands_PTV_cap100_p01_m001.rds",
             "Combined_Probands_PTV_cap100_p01_m003.rds","Proband_PTV_cap250_p01_m001.rds",
             "Proband_PTV_cap250_p01_m003.rds")
  names(paths) <- c("Combined Probands PTV","Combined Probands Mis2","Developmental disorders PTV",
                    "Developmental disorders Mis2")
  ces <- unique(read.delim(file.path(root,"inputs/gene_exclusions/updated_CES_genes_TableS36.tsv"))[[1]])
  map <- read.delim(file.path(root,"inputs/constraint/gnomad.v4.1.constraint_metrics.tsv"))
  stopifnot(all(ces %in% map$gene))
  ids <- unique(map$gene_id[map$gene %in% ces])
  xs <- list(); reports <- list()
  for(i in seq_along(paths)) {
    a <- readRDS(file.path(root,file.path("jobs", Sys.getenv("RUN_DATE"), "main"),paths[i]))$arguments
    keep <- !rownames(a$input_data) %in% ids
    d <- a$input_data[keep,,drop=FALSE]; f <- a$features[keep,,drop=FALSE]
    set.seed(if(i>=3) 20261104L else 20261103L)
    draws <- replicate(100,sample.int(nrow(d),nrow(d),replace=TRUE))
    xs[[names(paths)[i]]] <- list(input=d,features=f,draws=draws,
      prevalence=a$prevalence,scale=c(asd_scale,1,dd_scale,1)[i])
    reports[[i]] <- data.frame(model=names(paths)[i],gene_id=rownames(a$input_data)[!keep])
  }
  validate_upper_bound_pairs(xs)
  caps <- sort(unique(c(1/xs[[1]]$prevalence,60,80,100,125,150,200,250,300,400)))
  inventory <- expand.grid(RR_cap=caps,model=names(xs),stringsAsFactors=FALSE)
  saveRDS(list(models=xs,inventory=inventory),workspace)
  write.table(inventory,file.path(out,"inventory.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  write.table(do.call(rbind,reports),file.path(out,"excluded_genes.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
} else if(stage=="worker") {
  w <- readRDS(workspace); id <- as.integer(args[2]); stopifnot(id %in% seq_len(nrow(w$inventory)))
  job <- w$inventory[id,]; x <- w$models[[job$model]]
  path <- file.path(out,sprintf("job_%02d.rds",id))
  if(file.exists(path)) stop("Refusing to overwrite completed sensitivity job")
  # Match the primary analyses: ten midpoint integration points per component.
  t <- initialize_model(poisson_uniform_likelihood,x$input,seq(0,log(job$RR_cap),length.out=10),x$features,10)
  q <- moments(t,x$input); rows <- list(); ll <- NA_real_
  for(i in 0:100) {
    idx <- if(i==0) seq_len(nrow(x$input)) else x$draws[,i]
    m <- subset_model_likelihood_rows(t,idx); m$features <- t$features[idx,,drop=FALSE]
    m <- MixSQP_fit(m,allow_status_return=TRUE)
    stopifnot(isTRUE(m$fit_status$converged),isTRUE(m$fit_status$usable))
    if(i==0) ll <- absolute_mixture_log_likelihood(m)
    v <- measure(m,x$input[idx,,drop=FALSE],subset_moments(q,idx),x$prevalence,x$scale)
    stopifnot(all(is.finite(v)))
    rows[[i+1]] <- data.frame(model=job$model,RR_cap=job$RR_cap,replicate=i,metric=names(v),estimate=unname(v))
  }
  atomic_job_rds(list(job=job,draws=do.call(rbind,rows),log_likelihood=ll),path)
} else if(stage=="collect") {
  w <- readRDS(workspace)
  z <- lapply(seq_len(nrow(w$inventory)),function(i) {
    a <- readRDS(file.path(out,sprintf("job_%02d.rds",i)))
    stopifnot(identical(a$job$model,w$inventory$model[i]),
              identical(a$job$RR_cap,w$inventory$RR_cap[i])); a
  })
  d <- do.call(rbind,lapply(z,`[[`,"draws"))
  validate_upper_bound_pairs(w$models)
  stopifnot(nrow(d)==nrow(w$inventory)*101L*3L,!anyDuplicated(d[c("model","RR_cap","replicate","metric")]))
  write.table(d,file.path(out,"full_metric_class_draws.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  points <- do.call(rbind,lapply(z,function(a) data.frame(a$job,log_likelihood=a$log_likelihood)))
  points$dataset <- ifelse(grepl("^Developmental disorders ",points$model),"DD: PTV + Mis2","ASD: PTV + Mis2")
  pooled <- aggregate(log_likelihood~RR_cap+dataset,points,sum)
  pooled$loglik_difference_from_best <- ave(pooled$log_likelihood,pooled$dataset,FUN=function(x)x-max(x))
  write.table(pooled,file.path(out,"pooled_estimates.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
all<-pool_upper_bound_draws(d)
s<-subset(all,replicate==0);s$replicate<-NULL
for(i in seq_len(nrow(s))) {
  v<-subset(all,replicate>0&RR_cap==s$RR_cap[i]&metric==s$metric[i]&dataset==s$dataset[i])$estimate
  stopifnot(length(v)==100)
  s$lower95[i]<-quantile(v,.025);s$upper95[i]<-quantile(v,.975)
}
ll<-read.delim(file.path(out,"pooled_estimates.tsv"))
ll<-data.frame(RR_cap=ll$RR_cap,metric="loglik",estimate=ll$loglik_difference_from_best,
  dataset=ifelse(grepl("ASD",ll$dataset),"ASD","DD"),lower95=NA_real_,upper95=NA_real_)
s<-rbind(s,ll[names(s)])
write.table(s,file.path(out,"full_metric_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
library(ggplot2)
labels<-c(loglik="Log likelihood relative to best",mutvar="Mutational variance (%)",
  fraction_cases_RR5="Cases with RR > 5 variants (%)",effective_penetrance="PTV effective penetrance (%)")
s$value<-s$estimate*ifelse(s$metric=="loglik",1,100)
s$metric<-factor(s$metric,levels=names(labels),labels=labels)
s$primary_cap<-ifelse(s$dataset=="ASD",100,250)
p<-ggplot(s,aes(RR_cap,value,color=dataset))+
  geom_ribbon(data=subset(s,!is.na(lower95)),aes(ymin=100*lower95,ymax=100*upper95,fill=dataset),alpha=.18,color=NA)+
  geom_vline(aes(xintercept=primary_cap),color="grey65",linetype="dashed")+
  geom_line(linewidth=.6)+geom_point(size=1.8)+
  facet_wrap(vars(metric,dataset),scales="free_y",ncol=2)+
  scale_color_manual(values=c("#375E97","#B66C32"))+
  scale_fill_manual(values=c("#375E97","#B66C32"))+
  scale_x_continuous(breaks=c(45,100,200,300,400))+
  labs(x="Upper rate-ratio bound",y=NULL)+theme_bw(base_size=11,base_family="Arial")+
  theme(legend.position="none",panel.grid.minor=element_blank(),strip.background=element_rect(fill="grey95"))
ggsave(file.path(out,"UpperRRSweep_all_metrics.pdf"),p,width=10,height=11,device=cairo_pdf)
ggsave(file.path(out,"UpperRRSweep_all_metrics.png"),p,width=10,height=11,dpi=180)
cat(sprintf("Finished: %d point fits and %d bootstrap results, all three metrics.\n",
            nrow(w$inventory),100L*nrow(w$inventory)))

} else stop("Expected prepare, worker or collect")
