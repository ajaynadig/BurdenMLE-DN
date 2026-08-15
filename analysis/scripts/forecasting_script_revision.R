library(SummarizedExperiment)

# Forecast gene discovery from the fitted PTV effect-size distributions.
# Plotting is intentionally separate so the established manuscript ggplot code
# can be retained without mixing it with simulation logic.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit)) return(sub(paste0("^", flag, "="), "", hit[1]))
  i <- match(flag, args)
  if (!is.na(i)) return(args[i + 1])
  default
}

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_file <- normalizePath(script_file, mustWork = TRUE)
repo_dir <- dirname(dirname(script_file))
final_runs_dir <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT", unset = file.path(repo_dir, "analysis")), mustWork = FALSE)

model_candidates <- sort(list.files(
  file.path(final_runs_dir, "outputs", "data"),
  pattern = "^models_autism_mixsqp_.*\\.Rdata$", full.names = TRUE
))
model_file <- get_arg(
  "--model-file",
  if (length(model_candidates)) tail(model_candidates, 1) else NA_character_
)
output_dir <- get_arg(
  "--output-dir", file.path(final_runs_dir, "outputs", "derived", "forecasting")
)
num_iter <- as.integer(get_arg("--iterations", 100))
n_points <- as.integer(get_arg("--points", 20))
seed <- as.integer(get_arg("--seed", 24312342))
cores <- as.integer(get_arg("--cores", 1))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (is.na(model_file) || !file.exists(model_file)) stop("Main autism model not found.")

source_files <- c(
  "BurdenMLE_DN.R", "estimate_mutvar.R", "io.R", "likelihoods.R",
  "model.R", "secondary_analysis_functions.R"
)
for (source_file in source_files) source(if (identical(source_file, "secondary_analysis_functions.R")) file.path(repo_dir, "analysis", "scripts", source_file) else file.path(repo_dir, "R", source_file))
load(model_file)

model_index <- setNames(seq_along(autism_data$loop_vars$names), autism_data$loop_vars$names)
required <- c(
  "Combined Probands PTV", "Combined SPARK PTV", "Combined ASC PTV",
  "Combined GeneDx PTV", "Combined DDID PTV", "Combined Non-DDID PTV"
)
if (anyNA(model_index[required])) {
  stop("Missing forecasting models: ", paste(required[is.na(model_index[required])], collapse = ", "))
}
models <- BurdenMLE_DN_models_autism[[1]]
get_model <- function(name) models[[model_index[[name]]]]
get_data <- function(name) get_genetic_data(model_index[[name]], autism_data)$genetic_data

overall_name <- "Combined Probands PTV"
overall_model <- get_model(overall_name)
overall_data <- get_data(overall_name)

# Cohort-specific effect scaling uses the overall posterior distribution and
# directly observed cohort counts, matching the manuscript methods.
cohort_names <- c("SPARK", "ASC", "GeneDx")
cohort_scaling <- setNames(vapply(cohort_names, function(cohort) {
  binomial_analysis(
    model = overall_model,
    genetic_data_total = overall_data,
    genetic_data_subsample = get_data(paste("Combined", cohort, "PTV"))
  )$MLE
}, numeric(1)), cohort_names)

forecast_once <- function(old_data, sampler, n_new, gamma_scaling_factor,
                          thresholds = c(0, 2, 5, 10, 20)) {
  log_rr <- posterior_gene_samples(sampler, num_samples = 1L)[, 1L]
  new_counts <- rpois(
    nrow(old_data),
    2 * n_new * exp(log_rr * gamma_scaling_factor) * old_data$case_rate
  )
  combined_counts <- old_data$case_count + new_counts
  total_n <- old_data$N[1] + n_new
  expected_count <- 2 * total_n * old_data$case_rate
  poisson_p <- ppois(combined_counts - 0.001, expected_count, lower.tail = FALSE)
  significant_fdr <- p.adjust(poisson_p, method = "fdr") < 0.05
  significant_bonferroni <- poisson_p < 0.05 / length(poisson_p)
  true_rr <- exp(log_rr)

  data.frame(
    threshold = thresholds,
    count_fdr = vapply(thresholds, function(x) sum(significant_fdr & true_rr >= x), numeric(1)),
    count_bonferroni = vapply(
      thresholds, function(x) sum(significant_bonferroni & true_rr >= x), numeric(1)
    )
  )
}

overall_sampler <- posterior_gene_sampler(overall_model, overall_data)
ddid_sampler <- posterior_gene_sampler(
  get_model("Combined DDID PTV"), get_data("Combined DDID PTV")
)
non_ddid_sampler <- posterior_gene_sampler(
  get_model("Combined Non-DDID PTV"), get_data("Combined Non-DDID PTV")
)

scenario_list <- list(
  SPARK = list(data = overall_data, sampler = overall_sampler,
               scaling = cohort_scaling[["SPARK"]], group = "cohort"),
  ASC = list(data = overall_data, sampler = overall_sampler,
             scaling = cohort_scaling[["ASC"]], group = "cohort"),
  GeneDx = list(data = overall_data, sampler = overall_sampler,
                scaling = cohort_scaling[["GeneDx"]], group = "cohort"),
  # Both phenotype scenarios begin with the same complete observed dataset.
  # Only the distribution generating newly added cases differs.
  DDID = list(data = overall_data,
              sampler = ddid_sampler,
              scaling = 1, group = "ddid"),
  `Non-DDID` = list(data = overall_data,
                    sampler = non_ddid_sampler,
                    scaling = 1, group = "ddid")
)

if (!identical(rownames(overall_data), rownames(get_data("Combined DDID PTV"))) ||
    !identical(rownames(overall_data), rownames(get_data("Combined Non-DDID PTV")))) {
  stop("Overall and DDID-stratified gene orders differ; align genes before forecasting.")
}

max_new_cases <- as.integer(get_arg("--max-new-cases", overall_data$N[1]))
min_new_cases <- as.integer(get_arg("--min-new-cases", 1000))
n_positive_points <- max(1L, n_points - 1L)
n_new_range <- unique(c(
  0L,
  round(seq(min_new_cases, max_new_cases, length.out = n_positive_points))
))
tasks <- expand.grid(
  n_new_case = n_new_range,
  iteration = seq_len(num_iter),
  scenario = names(scenario_list),
  stringsAsFactors = FALSE
)

run_task <- function(task_index) {
  task <- tasks[task_index, ]
  scenario <- scenario_list[[task$scenario]]
  # A task-specific seed makes results invariant to core count and task order.
  set.seed(seed + task_index * 104729L)
  result <- forecast_once(
    scenario$data, scenario$sampler, task$n_new_case, scenario$scaling
  )
  result$scenario <- task$scenario
  result$scenario_group <- scenario$group
  result$n_new_case <- task$n_new_case
  result$total_sample_size <- scenario$data$N[1] + task$n_new_case
  result$iteration <- task$iteration
  result
}

cat("Forecasting", nrow(tasks), "scenario/iteration/sample-size combinations...\n")
if (cores > 1 && .Platform$OS.type != "windows") {
  task_results <- parallel::mclapply(seq_len(nrow(tasks)), run_task, mc.cores = cores)
} else {
  progress <- txtProgressBar(min = 0, max = nrow(tasks), style = 3)
  task_results <- vector("list", nrow(tasks))
  for (i in seq_len(nrow(tasks))) {
    task_results[[i]] <- run_task(i)
    setTxtProgressBar(progress, i)
  }
  close(progress)
}
forecast_full <- do.call(rbind, task_results)

summarize_forecast <- function(data) {
  keys <- unique(data[c("threshold", "scenario", "n_new_case", "total_sample_size")])
  out <- lapply(seq_len(nrow(keys)), function(i) {
    rows <- data$threshold == keys$threshold[i] &
      data$scenario == keys$scenario[i] &
      data$n_new_case == keys$n_new_case[i]
    fdr <- data$count_fdr[rows]
    bonf <- data$count_bonferroni[rows]
    cbind(
      keys[i, ],
      data.frame(
        count_fdr_mean = mean(fdr), count_fdr_sd = sd(fdr),
        count_fdr_lower = unname(quantile(fdr, 0.025)),
        count_fdr_upper = unname(quantile(fdr, 0.975)),
        count_bonferroni_mean = mean(bonf), count_bonferroni_sd = sd(bonf),
        count_bonferroni_lower = unname(quantile(bonf, 0.025)),
        count_bonferroni_upper = unname(quantile(bonf, 0.975))
      )
    )
  })
  do.call(rbind, out)
}

cohort_full <- forecast_full[forecast_full$scenario_group == "cohort", ]
ddid_full <- forecast_full[forecast_full$scenario_group == "ddid", ]
cohort_summary <- summarize_forecast(cohort_full)
ddid_summary <- summarize_forecast(ddid_full)

write.table(cohort_full, file.path(output_dir, "cohort_forecast_full.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cohort_summary, file.path(output_dir, "cohort_forecast_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ddid_full, file.path(output_dir, "ddid_forecast_full.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ddid_summary, file.path(output_dir, "ddid_forecast_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(
  list(
    metadata = list(model_file = normalizePath(model_file), seed = seed,
                    iterations = num_iter, points = length(n_new_range),
                    cohort_scaling = cohort_scaling),
    cohort = cohort_summary, ddid = ddid_summary
  ),
  file.path(output_dir, "forecast_summary.rds")
)

cat("Cohort effect-size scaling factors:\n")
print(cohort_scaling)
cat("Forecast outputs written to", output_dir, "\n")
