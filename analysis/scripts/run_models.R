# Command-line usage:
# Rscript run_models.R --mode test
# Rscript run_models.R --mode main --run-autism true --run-ddd false
library(pbapply)
library(data.table)
library(ggplot2)
library(patchwork)
library(SummarizedExperiment)

args <- commandArgs(trailingOnly = TRUE)

script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(script_file)) {
  script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
}
script_file <- normalizePath(script_file, mustWork = TRUE)
repo_dir <- dirname(dirname(dirname(script_file)))
source(file.path(repo_dir, "analysis", "scripts", "run_models_cli.R"))
source(file.path(repo_dir, "analysis", "scripts", "model_artifacts.R"))
source(file.path(repo_dir, "analysis", "scripts", "model_jobs.R"))
job_stage <- Sys.getenv("BURDENMLEDN_JOB_STAGE", "serial")
if (!job_stage %in% c("serial", "prepare", "collect")) stop("Invalid BURDENMLEDN_JOB_STAGE")
config <- parse_run_models_args(args)
if (isTRUE(config$help)) {
  cat(run_models_help())
  quit(status = 0)
}
analysis_mode <- config$analysis_mode
run_autism <- config$run_autism
run_ddd <- config$run_ddd
n_bootstraps <- config$n_bootstraps
random_seed <- config$random_seed
run_date <- config$run_date
autism_max_effect_size <- config$autism_max_effect_size
ddd_max_effect_size <- config$ddd_max_effect_size
ces_gene_file <- config$ces_gene_file
optimizer <- config$optimizer
set.seed(random_seed)

final_runs_dir <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT", unset = file.path(repo_dir, "analysis")), mustWork = FALSE)
input_dir <- file.path(final_runs_dir, "inputs")
output_data_dir <- file.path(final_runs_dir, "outputs", "data")
output_figure_dir <- file.path(final_runs_dir, "outputs", "figures", "diagnostics", "model_tests")
dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_figure_dir, recursive = TRUE, showWarnings = FALSE)
job_dir <- file.path(final_runs_dir, "jobs", run_date, analysis_mode)
job_paths <- character()
bootstrap_group <- 0L
# Preparing writes only job specifications, never placeholder fit artifacts.
save_model_artifact <- function(..., file) {
  if (job_stage == "prepare") return(invisible(NULL))
  if (file.exists(file)) stop("Refusing to overwrite model artifact: ", file)
  objects <- as.list(substitute(list(...)))[-1L]
  object_names <- vapply(objects, deparse, character(1))
  tmp <- tempfile(".partial-", tmpdir = dirname(file))
  on.exit(unlink(tmp), add = TRUE)
  save(list = object_names, file = tmp, envir = parent.frame())
  if (!file.rename(tmp, file)) stop("Cannot publish model artifact: ", file)
}

source_files <- c(
  "BurdenMLE_DN.R", "EM.R", "MixSQP.R", "estimate_mutvar.R", "io.R",
  "likelihoods.R", "model.R", "secondary_analysis_functions.R"
)
for (source_file in source_files) {
  source(if (identical(source_file, "secondary_analysis_functions.R")) file.path(repo_dir, "analysis", "scripts", source_file) else file.path(repo_dir, "R", source_file))
}

fit_BurdenMLE_DN_models <- function(data,
                                model_names,
                                prev_factor_indices,
                                max_effect_size,
                                bootstrap_count,
                                optimizer,
                                genes_to_keep = rownames(data$counts),
                                bootstrap_source = NULL) {
  missing_models <- setdiff(model_names, data$loop_vars$names)
  if (length(missing_models) > 0) {
    stop("Requested models are absent from the setup: ",
         paste(missing_models, collapse = ", "))
  }
  if (!all(prev_factor_indices %in% seq_along(data$prev_factors))) {
    stop("Invalid prevalence-factor index requested.")
  }

  model_indices <- match(model_names, data$loop_vars$names)
  bootstrap_gene_ids <- rownames(data$counts)[
    rownames(data$counts) %in% genes_to_keep
  ]
  bootstrap_group <<- bootstrap_group + 1L
  group_seed <- as.double(random_seed) +
    1000L * match(analysis_mode, c("main", "no_ces", "no_overlap", "age", "test")) + bootstrap_group
  if (group_seed > .Machine$integer.max) stop("Bootstrap group seed exceeds integer range.")
  bootstrap_samples <- if (job_stage == "collect") NULL else
    model_job_bootstraps(bootstrap_gene_ids, bootstrap_count, as.integer(group_seed))
  saved_draw_dir <- Sys.getenv("BURDENMLEDN_PREVIOUS_BOOTSTRAPS", "")
  if (nzchar(saved_draw_dir) && job_stage != "collect") {
    if (is.null(bootstrap_source)) stop("No previous-draw source assigned to this model group.")
    saved <- readRDS(file.path(saved_draw_dir, paste0(bootstrap_source, ".rds")))
    if (!identical(saved$gene_ids, bootstrap_gene_ids) ||
        !identical(dim(saved$bootstrap_samples), c(length(bootstrap_gene_ids), as.integer(bootstrap_count))) ||
        anyNA(saved$bootstrap_samples) || !all(saved$bootstrap_samples %in% bootstrap_gene_ids)) {
      stop("Previous bootstrap draws do not match this model's ordered gene universe.")
    }
    bootstrap_samples <- saved$bootstrap_samples
    message("Using exact saved bootstrap gene draws: ", bootstrap_source)
  }

  fit_model_flags <- if (is.null(data$loop_vars$fit_model)) {
    rep(TRUE, length(data$loop_vars$names))
  } else {
    data$loop_vars$fit_model
  }
  eligible_model_indices <- intersect(model_indices, which(fit_model_flags))
  total_fits <- length(prev_factor_indices) * length(eligible_model_indices)
  completed_fits <- 0L
  progress_bar <- txtProgressBar(min = 0, max = total_fits, initial = 0, style = 3)
  on.exit(close(progress_bar), add = TRUE)

  lapply(seq_along(data$prev_factors), function(prev_factor_index) {
    prev_factor <- data$prev_factors[prev_factor_index]
    print(paste("Prevalence scaling factor:", prev_factor))

    lapply(seq_along(data$loop_vars$names), function(i) {
      should_fit <- prev_factor_index %in% prev_factor_indices &&
        i %in% model_indices &&
        (is.null(data$loop_vars$fit_model) || data$loop_vars$fit_model[i])

      if (!should_fit) {
        return(NA)
      }

      cat(
        "\n",
        sprintf("Model %d/%d: %s", completed_fits + 1L, total_fits,
                data$loop_vars$names[i]),
        "\n",
        sep = ""
      )
      processed_input <- get_genetic_data(i, data)
      gene_keep <- rownames(processed_input$genetic_data) %in% genes_to_keep
      input_df <- processed_input$genetic_data[gene_keep, ]
      features <- processed_input$features[gene_keep, ]
      print(colSums(features))

      fit_arguments <- list(
        input_data = input_df,
        component_endpoints = seq(
          0, log(max_effect_size), length.out = 10
        ),
        features = features,
        null_sim = FALSE,
        prevalence = data$loop_vars$prevalences[i] * prev_factor,
        estimate_posteriors = TRUE,
        bootstrap_samples = bootstrap_samples,
        n_boot = bootstrap_count,
        optimizer = optimizer
      )

      if (job_stage == "serial") {
        model <- do.call(BurdenMLE_DN, fit_arguments)
      } else {
        collection <- gsub("[^A-Za-z0-9]+", "_", data$loop_vars$names[1])
        job_path <- file.path(job_dir, sprintf("%s_cap%s_p%02d_m%03d.rds",
                                               collection, max_effect_size, prev_factor_index, i))
        if (job_stage == "collect") {
          # Reuse the actual draws that produced the checkpoint. Never redraw here.
          fit_arguments$bootstrap_samples <- readRDS(job_path)$arguments$bootstrap_samples
          if (!identical(dim(fit_arguments$bootstrap_samples),
                         c(length(bootstrap_gene_ids), as.integer(bootstrap_count))) ||
              !all(fit_arguments$bootstrap_samples %in% bootstrap_gene_ids)) {
            stop("Invalid saved bootstrap samples: ", job_path)
          }
        }
        prepare_model_job(fit_arguments, data$loop_vars$names[i], job_path, repo_dir)
        job_paths <<- c(job_paths, job_path)
        if (job_stage == "prepare") return(NA)
        model <- read_model_job_result(job_path, repo_dir)
      }

      print(model$mutvar_output$total_mutvar)
      print(model$mutvar_output$enrichment)
      print(model$mutvar_output$mutvar_CI)
      completed_fits <<- completed_fits + 1L
      setTxtProgressBar(progress_bar, completed_fits)
      model
    })
  })
}

current_date <- run_date
optimizer_suffix <- if (optimizer == "EM") "" else paste0("_", tolower(optimizer))
artifact_paths <- character()

setup_path <- file.path(job_dir, "setup_snapshot.rds")
use_snapshot <- job_stage != "serial" && file.exists(setup_path)
if (job_stage == "collect" && !use_snapshot) stop("Collection requires a prepared setup snapshot.")
if (use_snapshot) {
  list2env(read_model_setup(setup_path, config, repo_dir), envir = .GlobalEnv)
} else {

if ((analysis_mode == "main" && run_autism) ||
    analysis_mode %in% c("no_ces", "test")) {
  source(file.path(repo_dir, "analysis", "scripts", "set_up_asc.R"))
}

if (analysis_mode == "no_overlap") {
  source(file.path(repo_dir, "analysis", "scripts", "set_up_asc_NoKaplanisOverlap.R"))
}

if (analysis_mode == "age") {
  source(file.path(repo_dir, "analysis", "scripts", "set_up_SPARK_agephenos.R"))
}

if (analysis_mode == "main" && run_ddd) {
  source(file.path(repo_dir, "analysis", "scripts", "set_up_kaplanis.R"))
}

if (job_stage == "prepare") {
  object_names <- intersect(c("autism_data", "kaplanis_data", "autism_counts",
                              "age_phenotypes_spark", "theme_bhr_legend_gridlines"), ls(envir = .GlobalEnv))
  atomic_job_rds(list(config = config, code = model_setup_code(repo_dir),
                     objects = mget(object_names, envir = .GlobalEnv)), setup_path)
  writeLines(unname(tools::md5sum(setup_path)), paste0(setup_path, ".md5"))
}
}

if (analysis_mode == "main") {
  if (run_autism) {
    autism_model_names <- autism_data$loop_vars$names[
      autism_data$loop_vars$fit_model
    ]
    BurdenMLE_DN_models_autism <- fit_BurdenMLE_DN_models(
      data = autism_data,
      model_names = autism_model_names,
      bootstrap_source = "main_autism",
      prev_factor_indices = seq_along(autism_data$prev_factors),
      max_effect_size = autism_max_effect_size,
      bootstrap_count = n_bootstraps,
      optimizer = optimizer
    )
    autism_artifact_path <- file.path(
      output_data_dir,
      paste0("models_autism", optimizer_suffix, "_", current_date, ".Rdata")
    )
    save_model_artifact(
      BurdenMLE_DN_models_autism, autism_data,
      file = autism_artifact_path
    )
    artifact_paths["autism"] <- autism_artifact_path
  }

  if (run_ddd) {
    BurdenMLE_DN_models_DDD <- fit_BurdenMLE_DN_models(
      data = kaplanis_data,
      model_names = kaplanis_data$loop_vars$names,
      bootstrap_source = "main_ddd",
      prev_factor_indices = seq_along(kaplanis_data$prev_factors),
      max_effect_size = ddd_max_effect_size,
      bootstrap_count = n_bootstraps,
      optimizer = optimizer
    )
    ddd_artifact_path <- file.path(
      output_data_dir,
      paste0("models_ddd", optimizer_suffix, "_", current_date, ".Rdata")
    )
    save_model_artifact(
      BurdenMLE_DN_models_DDD, kaplanis_data,
      file = ddd_artifact_path
    )
    artifact_paths["ddd"] <- ddd_artifact_path
  }
}

if (analysis_mode == "age") {
  age_model_names <- autism_data$loop_vars$names
  BurdenMLE_DN_models_autism_AgePhenos <- fit_BurdenMLE_DN_models(
    data = autism_data,
    model_names = age_model_names,
    bootstrap_source = "age",
    prev_factor_indices = 1,
    max_effect_size = autism_max_effect_size,
    bootstrap_count = n_bootstraps,
    optimizer = optimizer
  )
  autism_data_AgePhenos <- autism_data
  autism_counts_AgePhenos <- autism_counts
  age_artifact_path <- file.path(
    output_data_dir,
    paste0("models_autism_AgePhenos", optimizer_suffix, "_", current_date, ".Rdata")
  )
  save_model_artifact(
    BurdenMLE_DN_models_autism_AgePhenos,
    autism_data_AgePhenos,
    autism_counts_AgePhenos,
    age_phenotypes_spark,
    file = age_artifact_path
  )
  artifact_paths["age"] <- age_artifact_path
}

if (analysis_mode == "no_ces") {
  no_ces_model_names <- c(
    "Combined Probands PTV", "Combined Siblings PTV",
    "Combined Probands Mis2", "Combined Siblings Mis2",
    "Combined Probands Mis1", "Combined Siblings Mis1",
    "Combined Probands Mis0", "Combined Siblings Mis0",
    "Combined Probands Syn", "Combined Siblings Syn"
  )

  gnomad_v2 <- data.frame(fread(
    file.path(input_dir, "constraint", "gnomad.v2.1.1.lof_metrics.by_gene.txt")
  ))
  if (is.null(ces_gene_file)) {
    ces_gene_file <- file.path(
      input_dir, "gene_exclusions", "updated_CES_genes_TableS36.tsv"
    )
  }
  ces_gene_file <- normalizePath(ces_gene_file, mustWork = TRUE)
  ces_gene_table <- read.delim(
    ces_gene_file, header = TRUE, stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (ncol(ces_gene_table) < 1L) stop("CES gene table has no columns.")
  ces_gene_symbols <- trimws(as.character(ces_gene_table[[1L]]))
  if (anyNA(ces_gene_symbols) || any(!nzchar(ces_gene_symbols)) ||
      anyDuplicated(ces_gene_symbols)) {
    stop("The first CES gene-table column must contain unique, nonempty symbols.")
  }
  lookup_symbols <- ces_gene_symbols
  lookup_symbols[lookup_symbols == "CERT1"] <- "COL4A3BP"
  lookup_symbols[lookup_symbols == "NARS1"] <- "NARS"
  ces_gene_ids <- gnomad_v2$gene_id[
    match(lookup_symbols, gnomad_v2$gene)
  ]
  if (anyNA(ces_gene_ids)) {
    stop("CES symbols absent from the gene reference: ",
         paste(ces_gene_symbols[is.na(ces_gene_ids)], collapse = ", "))
  }
  genes_to_keep <- setdiff(rownames(autism_data$counts), ces_gene_ids)
  print(paste("Removing", sum(!rownames(autism_data$counts) %in% genes_to_keep),
              "CES genes"))

  BurdenMLE_DN_models_autism_noCES <- fit_BurdenMLE_DN_models(
    data = autism_data,
    model_names = no_ces_model_names,
    bootstrap_source = "no_ces",
    prev_factor_indices = 1,
    max_effect_size = autism_max_effect_size,
    bootstrap_count = n_bootstraps,
    optimizer = optimizer,
    genes_to_keep = genes_to_keep
  )
  no_ces_artifact_path <- file.path(
    output_data_dir,
    paste0("models_autism_noCES", optimizer_suffix, "_", current_date, ".Rdata")
  )
  save_model_artifact(
    BurdenMLE_DN_models_autism_noCES, autism_data, ces_gene_ids,
    file = no_ces_artifact_path
  )
  artifact_paths["no_ces"] <- no_ces_artifact_path
}

if (analysis_mode == "no_overlap") {
  no_overlap_model_names <- autism_data$loop_vars$names
  BurdenMLE_DN_models_autism_NoKaplanis <- fit_BurdenMLE_DN_models(
    data = autism_data,
    model_names = no_overlap_model_names,
    bootstrap_source = "no_overlap",
    prev_factor_indices = seq_along(autism_data$prev_factors),
    max_effect_size = autism_max_effect_size,
    bootstrap_count = n_bootstraps,
    optimizer = optimizer
  )
  autism_data_NoKaplanis <- autism_data
  autism_counts_NoKaplanis <- autism_counts
  no_overlap_artifact_path <- file.path(
    output_data_dir,
    paste0("models_autism_no_overlap", optimizer_suffix, "_", current_date, ".Rdata")
  )
  save_model_artifact(
    BurdenMLE_DN_models_autism_NoKaplanis,
    autism_data_NoKaplanis,
    autism_counts_NoKaplanis,
    file = no_overlap_artifact_path
  )
  artifact_paths["no_overlap"] <- no_overlap_artifact_path
}

if (analysis_mode == "test") {
  test_model_names <- c(
    "Combined Probands PTV",
    "Combined Probands Mis2",
    "Combined Probands Mis1",
    "Combined Probands Mis0"
  )
  BurdenMLE_DN_models_autism_test <- fit_BurdenMLE_DN_models(
    data = autism_data,
    model_names = test_model_names,
    prev_factor_indices = 1,
    max_effect_size = autism_max_effect_size,
    bootstrap_count = n_bootstraps,
    optimizer = optimizer
  )
  test_artifact_path <- file.path(
    output_data_dir,
    paste0("models_autism_test", optimizer_suffix, "_", current_date, ".Rdata")
  )
  save_model_artifact(
    BurdenMLE_DN_models_autism_test, autism_data,
    file = test_artifact_path
  )
  artifact_paths["test"] <- test_artifact_path

  if (job_stage != "prepare") {
  test_summary <- mutvar_enrichment_table(
    autism_data, BurdenMLE_DN_models_autism_test
  )
  test_mutvar <- test_summary$mutvar[
    test_summary$mutvar$name %in% test_model_names &
      test_summary$mutvar$prev_factor == autism_data$prev_factors[1],
  ]
  test_enrichment <- test_summary$enrichment[
    test_summary$enrichment$name == "Combined Probands PTV" &
      test_summary$enrichment$prev_factor == autism_data$prev_factors[1] &
      test_summary$enrichment$annot != "LOEUF5",
  ]
  test_enrichment$annotation <- gsub("_", "\n", test_enrichment$annot)

  palette_variantclass <- c(
    "PTV" = "#FF5F8A", "Mis2" = "darkorange",
    "Mis1" = "#F9B332", "Mis0" = "#EFD09F"
  )
  test_mutvar_plot <- ggplot(
    test_mutvar,
    aes(x = factor(variant_class, levels = c("Mis0", "Mis1", "Mis2", "PTV")),
        y = mutvar, ymin = mutvar_lower, ymax = mutvar_upper, fill = variant_class)
  ) +
    geom_hline(yintercept = 0) +
    geom_pointrange(shape = 21, color = "black", size = 1) +
    scale_fill_manual(values = palette_variantclass) +
    theme_bhr_legend_gridlines() +
    labs(x = "Variant class", y = "Mutational variance\nObserved scale") +
    guides(fill = "none")

  test_enrichment_plot <- ggplot(
    test_enrichment,
    aes(x = factor(annotation, levels = annotation),
        y = mutvar_enrich, ymin = mutvar_enrich_lower, ymax = mutvar_enrich_upper)
  ) +
    geom_col(fill = "#FF5F8A", color = "black") +
    geom_errorbar(width = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme_bhr_legend_gridlines() +
    labs(x = "LOEUF/mutation-rate annotation", y = "PTV enrichment")

  test_figure <- test_mutvar_plot + test_enrichment_plot +
    plot_annotation(tag_levels = "A")
  ggsave(
    file.path(
      output_figure_dir,
      paste0("models_autism_test_diagnostic", optimizer_suffix, "_", current_date, ".pdf")
    ),
    test_figure, device = cairo_pdf, width = 11, height = 4.5
  )
  }
}

if (job_stage == "prepare") {
  writeLines(job_paths, file.path(job_dir, "tasks.txt"))
  cat("Prepared", length(job_paths), "model jobs in", job_dir, "\n")
  quit(status = 0)
}

manifest_path <- write_model_manifest(
  path = file.path(
    output_data_dir,
    paste0("model_manifest_", analysis_mode, "_", current_date, ".rds")
  ),
  run_id = current_date,
  mode = analysis_mode,
  repo_dir = repo_dir,
  settings = list(
    optimizer = optimizer,
    seed = random_seed,
    grid_size = 10L,
    bootstrap_count = n_bootstraps,
    max_effect_size = list(
      autism = autism_max_effect_size,
      ddd = ddd_max_effect_size
    )
  ),
  artifacts = artifact_paths
)
cat("Model manifest written to", manifest_path, "\n")
