library(SummarizedExperiment)

# Lightweight summaries for supplementary autism analyses whose estimates are
# already stored in the fitted model objects.

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
  "--model-file", if (length(model_candidates)) tail(model_candidates, 1) else NA_character_
)
output_dir <- get_arg(
  "--output-dir", file.path(final_runs_dir, "outputs", "derived", "main_autism")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (is.na(model_file) || !file.exists(model_file)) stop("Main autism model not found.")
for (source_file in c("io.R", "model.R", "secondary_analysis_functions.R")) {
  source(if (identical(source_file, "secondary_analysis_functions.R")) file.path(repo_dir, "analysis", "scripts", source_file) else file.path(repo_dir, "R", source_file))
}
load(model_file)
model_index <- setNames(seq_along(autism_data$loop_vars$names), autism_data$loop_vars$names)
variant_classes <- c("PTV", "Mis2", "Mis1", "Mis0", "Syn")

extract_mutvar <- function(prevalence_index, role) {
  model_names <- paste("Combined", role, variant_classes)
  if (anyNA(model_index[model_names])) {
    stop("Missing models: ", paste(model_names[is.na(model_index[model_names])], collapse = ", "))
  }
  do.call(rbind, lapply(seq_along(variant_classes), function(i) {
    model <- BurdenMLE_DN_models_autism[[prevalence_index]][[model_index[[model_names[i]]]]]
    scale <- if (variant_classes[i] == "PTV") autism_data$ptv_scale_factor else 1
    data.frame(
      role = role, variant_class = variant_classes[i],
      prevalence_factor = autism_data$prev_factors[prevalence_index],
      prevalence = autism_data$baseprev * autism_data$prev_factors[prevalence_index],
      estimate = model$mutvar_output$total_mutvar * scale,
      lower = model$mutvar_output$mutvar_CI[1] * scale,
      upper = model$mutvar_output$mutvar_CI[2] * scale
    )
  }))
}

sibling_mutvar <- extract_mutvar(1, "Siblings")
sibling_ptv_index <- model_index[["Combined Siblings PTV"]]
sibling_mutvar$n <- get_genetic_data(sibling_ptv_index, autism_data)$genetic_data$N[1]
prevalence_mutvar <- do.call(rbind, lapply(
  seq_along(autism_data$prev_factors), extract_mutvar, role = "Probands"
))

summary_output <- list(
  metadata = list(model_file = normalizePath(model_file), generated = format(Sys.time())),
  sibling_mutational_variance = sibling_mutvar,
  prevalence_mutational_variance = prevalence_mutvar
)
saveRDS(summary_output, file.path(output_dir, "secondary_autism_summary.rds"))
write.table(sibling_mutvar, file.path(output_dir, "sibling_mutational_variance.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(prevalence_mutvar, file.path(output_dir, "prevalence_mutational_variance.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Secondary autism summaries written to", output_dir, "\n")
print(sibling_mutvar, row.names = FALSE)
