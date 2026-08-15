# Command-line usage:
# Rscript run_models.R --mode test
# Rscript run_models.R --mode main --run-autism true --run-ddd false
library(pbapply)

args <- commandArgs(trailingOnly = TRUE)

script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(script_file)) {
  script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
}
script_file <- normalizePath(script_file, mustWork = TRUE)
repo_dir <- dirname(dirname(dirname(script_file)))
source(file.path(repo_dir, "analysis", "scripts", "run_models_cli.R"))
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
optimizer <- config$optimizer
set.seed(random_seed)

final_runs_dir <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT", unset = file.path(repo_dir, "analysis")), mustWork = FALSE)
input_dir <- file.path(final_runs_dir, "inputs")
output_data_dir <- file.path(final_runs_dir, "outputs", "data")
output_figure_dir <- file.path(final_runs_dir, "outputs", "figures", "diagnostics", "model_tests")
dir.create(output_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_figure_dir, recursive = TRUE, showWarnings = FALSE)

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
                                genes_to_keep = rownames(data$counts)) {
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
  bootstrap_samples <- matrix(
    vapply(seq_len(bootstrap_count), function(dummy) {
      bootstrap_gene_ids[
        sample(seq_along(bootstrap_gene_ids), replace = TRUE)
      ]
    }, character(length(bootstrap_gene_ids))),
    nrow = length(bootstrap_gene_ids),
    ncol = bootstrap_count
  )

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

      model <- BurdenMLE_DN(
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

if (analysis_mode == "main") {
  if (run_autism) {
    autism_model_names <- autism_data$loop_vars$names[
      autism_data$loop_vars$fit_model
    ]
    BurdenMLE_DN_models_autism <- fit_BurdenMLE_DN_models(
      data = autism_data,
      model_names = autism_model_names,
      prev_factor_indices = seq_along(autism_data$prev_factors),
      max_effect_size = autism_max_effect_size,
      bootstrap_count = n_bootstraps,
      optimizer = optimizer
    )
    save(
      BurdenMLE_DN_models_autism, autism_data,
      file = file.path(
        output_data_dir,
        paste0("models_autism", optimizer_suffix, "_", current_date, ".Rdata")
      )
    )
  }

  if (run_ddd) {
    BurdenMLE_DN_models_DDD <- fit_BurdenMLE_DN_models(
      data = kaplanis_data,
      model_names = kaplanis_data$loop_vars$names,
      prev_factor_indices = seq_along(kaplanis_data$prev_factors),
      max_effect_size = ddd_max_effect_size,
      bootstrap_count = n_bootstraps,
      optimizer = optimizer
    )
    save(
      BurdenMLE_DN_models_DDD, kaplanis_data,
      file = file.path(
        output_data_dir,
        paste0("models_ddd", optimizer_suffix, "_", current_date, ".Rdata")
      )
    )
  }
}

if (analysis_mode == "age") {
  age_model_names <- autism_data$loop_vars$names
  BurdenMLE_DN_models_autism_AgePhenos <- fit_BurdenMLE_DN_models(
    data = autism_data,
    model_names = age_model_names,
    prev_factor_indices = 1,
    max_effect_size = autism_max_effect_size,
    bootstrap_count = n_bootstraps,
    optimizer = optimizer
  )
  autism_data_AgePhenos <- autism_data
  autism_counts_AgePhenos <- autism_counts
  save(
    BurdenMLE_DN_models_autism_AgePhenos,
    autism_data_AgePhenos,
    autism_counts_AgePhenos,
    age_phenotypes_spark,
    file = file.path(
      output_data_dir,
      paste0("models_autism_AgePhenos", optimizer_suffix, "_", current_date, ".Rdata")
    )
  )
}

if (analysis_mode == "no_ces") {
  no_ces_model_names <- c(
    "Combined Probands PTV", "Combined Siblings PTV",
    "Combined Probands Mis2", "Combined Siblings Mis2",
    "Combined Probands Mis1", "Combined Siblings Mis1",
    "Combined Probands Mis0", "Combined Siblings Mis0",
    "Combined Probands Syn", "Combined Siblings Syn"
  )

  ces_gene_table <- read.csv(
    file.path(input_dir, "gene_exclusions", "Sepliyarskiy_SuppTable6.csv")
  )
  gnomad_v2 <- data.frame(fread(
    file.path(input_dir, "constraint", "gnomad.v2.1.1.lof_metrics.by_gene.txt")
  ))
  ces_gene_symbols <- ces_gene_table$GeneID[
    ces_gene_table$Category == "set 2"
  ]
  ces_gene_ids <- gnomad_v2$gene_id[
    match(ces_gene_symbols, gnomad_v2$gene)
  ]
  ces_gene_ids[ces_gene_symbols == "NARS1"] <- "ENSG00000134440"
  genes_to_keep <- setdiff(rownames(autism_data$counts), ces_gene_ids)
  print(paste("Removing", sum(!rownames(autism_data$counts) %in% genes_to_keep),
              "CES genes"))

  BurdenMLE_DN_models_autism_noCES <- fit_BurdenMLE_DN_models(
    data = autism_data,
    model_names = no_ces_model_names,
    prev_factor_indices = 1,
    max_effect_size = autism_max_effect_size,
    bootstrap_count = n_bootstraps,
    optimizer = optimizer,
    genes_to_keep = genes_to_keep
  )
  save(
    BurdenMLE_DN_models_autism_noCES, autism_data, ces_gene_ids,
    file = file.path(
      output_data_dir,
      paste0("models_autism_noCES", optimizer_suffix, "_", current_date, ".Rdata")
    )
  )
}

if (analysis_mode == "no_overlap") {
  no_overlap_model_names <- autism_data$loop_vars$names
  BurdenMLE_DN_models_autism_NoKaplanis <- fit_BurdenMLE_DN_models(
    data = autism_data,
    model_names = no_overlap_model_names,
    prev_factor_indices = seq_along(autism_data$prev_factors),
    max_effect_size = autism_max_effect_size,
    bootstrap_count = n_bootstraps,
    optimizer = optimizer
  )
  autism_data_NoKaplanis <- autism_data
  autism_counts_NoKaplanis <- autism_counts
  save(
    BurdenMLE_DN_models_autism_NoKaplanis,
    autism_data_NoKaplanis,
    autism_counts_NoKaplanis,
    file = file.path(
      output_data_dir,
      paste0("models_autism_no_overlap", optimizer_suffix, "_", current_date, ".Rdata")
    )
  )
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
  save(
    BurdenMLE_DN_models_autism_test, autism_data,
    file = file.path(
      output_data_dir,
      paste0("models_autism_test", optimizer_suffix, "_", current_date, ".Rdata")
    )
  )

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
