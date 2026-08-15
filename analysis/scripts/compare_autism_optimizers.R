library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default) {
  equals_match <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(equals_match) > 0) {
    return(sub(paste0("^", flag, "="), "", equals_match[1]))
  }
  flag_index <- match(flag, args)
  if (!is.na(flag_index)) {
    if (flag_index == length(args) || startsWith(args[flag_index + 1], "--")) {
      stop("Missing value after ", flag)
    }
    return(args[flag_index + 1])
  }
  default
}

script_file <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
)
script_file <- normalizePath(script_file, mustWork = TRUE)
repo_dir <- dirname(dirname(script_file))
final_runs_dir <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT", unset = file.path(repo_dir, "analysis")), mustWork = FALSE)
new_results_file <- get_arg(
  "--mixsqp-file",
  file.path(final_runs_dir, "outputs", "data",
            "models_autism_mixsqp_Jul16_26.Rdata")
)
old_results_file <- get_arg(
  "--em-file",
  file.path(final_runs_dir, "outputs", "data",
            "models_autism_em_Jul16_26.Rdata")
)
output_figure_dir <- file.path(final_runs_dir, "outputs", "figures", "diagnostics", "optimizer")
component_figure_dir <- file.path(output_figure_dir, "components")
output_table_dir <- file.path(final_runs_dir, "outputs", "tables", "diagnostics", "optimizer")
dir.create(output_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(component_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_table_dir, recursive = TRUE, showWarnings = FALSE)

for (path in c(new_results_file, old_results_file)) {
  if (!file.exists(path)) stop("Required file does not exist: ", path)
}

last_finite <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else tail(x, 1)
}

annotation_names <- function(model, n_annotations) {
  candidates <- c(
    names(model$mutvar_output$enrichment),
    rownames(model$delta),
    colnames(model$features)
  )
  candidates <- candidates[nzchar(candidates)]
  if (length(candidates) >= n_annotations) {
    return(candidates[seq_len(n_annotations)])
  }
  canonical <- c("LOEUF1_mu1", "LOEUF1_mu2", paste0("LOEUF", 2:5))
  if (n_annotations == length(canonical)) return(canonical)
  paste0("Annotation_", seq_len(n_annotations))
}

extract_results <- function(models, model_names, prevalence_factors,
                            method, ptv_scale_factor) {
  if (length(models) != length(prevalence_factors)) {
    stop(method, " prevalence dimension does not match its metadata.")
  }
  if (any(vapply(models, length, integer(1)) != length(model_names))) {
    stop(method, " model dimension does not match its reconstructed names.")
  }

  model_rows <- list()
  enrichment_rows <- list()
  row_number <- 0L
  enrichment_row_number <- 0L

  for (prevalence_index in seq_along(prevalence_factors)) {
    for (model_index in seq_along(model_names)) {
      model <- models[[prevalence_index]][[model_index]]
      if (!is.list(model) || is.null(model$mutvar_output)) next

      row_number <- row_number + 1L
      model_name <- model_names[model_index]
      variant_class <- if (grepl("PTV", model_name)) {
        "PTV"
      } else if (grepl("Mis2", model_name)) {
        "Mis2"
      } else if (grepl("Mis1", model_name)) {
        "Mis1"
      } else if (grepl("Mis0", model_name)) {
        "Mis0"
      } else if (grepl("Syn", model_name)) {
        "Syn"
      } else {
        "Other"
      }
      role <- if (grepl("Sibling", model_name)) "Sibling" else "Proband"
      mutvar_scale <- if (variant_class == "PTV") ptv_scale_factor else 1

      model_rows[[row_number]] <- data.frame(
        model_name = model_name,
        prevalence_index = prevalence_index,
        prevalence_factor = prevalence_factors[prevalence_index],
        variant_class = variant_class,
        role = role,
        optimizer = method,
        log_likelihood = model$fit_status$log_likelihood,
        mutational_variance =
          as.numeric(model$mutvar_output$total_mutvar) * mutvar_scale,
        mutational_variance_lower =
          as.numeric(model$mutvar_output$mutvar_CI[1]) * mutvar_scale,
        mutational_variance_upper =
          as.numeric(model$mutvar_output$mutvar_CI[2]) * mutvar_scale,
        effective_penetrance =
          as.numeric(model$penetrance$effective_penetrance),
        effective_penetrance_lower =
          as.numeric(model$penetrance$effective_penetrance_CI[1]),
        effective_penetrance_upper =
          as.numeric(model$penetrance$effective_penetrance_CI[2]),
        stringsAsFactors = FALSE
      )

      enrichment <- as.numeric(model$mutvar_output$enrichment)
      annotations <- annotation_names(model, length(enrichment))
      for (annotation_index in seq_along(enrichment)) {
        enrichment_row_number <- enrichment_row_number + 1L
        enrichment_rows[[enrichment_row_number]] <- data.frame(
          model_name = model_name,
          prevalence_index = prevalence_index,
          prevalence_factor = prevalence_factors[prevalence_index],
          variant_class = variant_class,
          role = role,
          annotation = annotations[annotation_index],
          optimizer = method,
          enrichment = enrichment[annotation_index],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  list(
    models = do.call(rbind, model_rows),
    enrichment = do.call(rbind, enrichment_rows)
  )
}

load_new_results <- function(path) {
  environment <- new.env(parent = emptyenv())
  loaded <- load(path, envir = environment)
  required <- c("BurdenMLE_DN_models_autism", "autism_data")
  if (!all(required %in% loaded)) {
    stop("MixSQP file must contain: ", paste(required, collapse = ", "))
  }
  data <- environment$autism_data
  extracted <- extract_results(
    environment$BurdenMLE_DN_models_autism,
    data$loop_vars$names,
    data$prev_factors,
    "MixSQP",
    data$ptv_scale_factor
  )
  metadata <- list(
    model_names = data$loop_vars$names,
    prevalence_factors = data$prev_factors,
    ptv_scale_factor = data$ptv_scale_factor
  )
  list(extracted = extracted, metadata = metadata)
}

cat("Loading and summarizing MixSQP results...\n")
new <- load_new_results(new_results_file)
gc()

cat("Loading and summarizing EM results...\n")
old_environment <- new.env(parent = emptyenv())
old_loaded <- load(old_results_file, envir = old_environment)
old_required <- c("BurdenMLE_DN_models_autism", "autism_data")
if (!all(old_required %in% old_loaded)) {
  stop("EM file must contain: ", paste(old_required, collapse = ", "))
}
old_models <- old_environment$BurdenMLE_DN_models_autism
old_data <- old_environment$autism_data
old <- extract_results(
  old_models,
  old_data$loop_vars$names,
  old_data$prev_factors,
  "EM",
  old_data$ptv_scale_factor
)
rm(old_environment, old_models, old_data)
gc()

model_keys <- c("model_name", "prevalence_index")
model_comparison <- merge(
  old$models,
  new$extracted$models,
  by = model_keys,
  suffixes = c("_em", "_mixsqp"),
  sort = FALSE
)
enrichment_comparison <- merge(
  old$enrichment,
  new$extracted$enrichment,
  by = c(model_keys, "annotation"),
  suffixes = c("_em", "_mixsqp"),
  sort = FALSE
)

if (nrow(model_comparison) == 0 || nrow(enrichment_comparison) == 0) {
  stop("No models matched between the EM and MixSQP results.")
}
if (any(model_comparison$prevalence_factor_em !=
        model_comparison$prevalence_factor_mixsqp)) {
  stop("Matched models have inconsistent prevalence factors.")
}

for (metric in c("log_likelihood", "mutational_variance",
                 "effective_penetrance")) {
  model_comparison[[paste0(metric, "_difference")]] <-
    model_comparison[[paste0(metric, "_mixsqp")]] -
    model_comparison[[paste0(metric, "_em")]]
}
for (metric in c("mutational_variance", "effective_penetrance")) {
  em_width <- model_comparison[[paste0(metric, "_upper_em")]] -
    model_comparison[[paste0(metric, "_lower_em")]]
  mixsqp_width <- model_comparison[[paste0(metric, "_upper_mixsqp")]] -
    model_comparison[[paste0(metric, "_lower_mixsqp")]]
  model_comparison[[paste0(metric, "_ci_width_em")]] <- em_width
  model_comparison[[paste0(metric, "_ci_width_mixsqp")]] <- mixsqp_width
  model_comparison[[paste0(metric, "_ci_width_ratio")]] <- mixsqp_width / em_width
}
enrichment_comparison$enrichment_difference <-
  enrichment_comparison$enrichment_mixsqp -
  enrichment_comparison$enrichment_em

write.table(
  model_comparison,
  file.path(output_table_dir, "autism_optimizer_comparison_models.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  enrichment_comparison,
  file.path(output_table_dir, "autism_optimizer_comparison_enrichment.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Keep the complete matched results in the audit tables, but focus the visual
# comparison on proband fits. Near-null sibling fits make normalized summaries
# visually unstable and are not central to the optimizer comparison.
model_plot_data <- model_comparison[
  model_comparison$role_em == "Proband" &
    model_comparison$prevalence_index == 1,
]
enrichment_plot_data <- enrichment_comparison[
  enrichment_comparison$role_em == "Proband" &
    enrichment_comparison$prevalence_index == 1,
]

theme_comparison <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      strip.background = element_rect(fill = "white")
    )
}

variant_colors <- c(
  PTV = "#D55E00", Mis2 = "#0072B2", Mis1 = "#009E73",
  Mis0 = "#E69F00", Syn = "#999999", Other = "#CC79A7"
)
annotation_colors <- c(
  LOEUF1_mu1 = "#440154FF", LOEUF1_mu2 = "#414487FF",
  LOEUF2 = "#2A788EFF", LOEUF3 = "#22A884FF",
  LOEUF4 = "#7AD151FF", LOEUF5 = "#FDE725FF"
)

comparison_plot <- function(data, x, y, title, axis_label,
                            color_variable = "variant_class_em",
                            color_values = variant_colors) {
  finite <- is.finite(data[[x]]) & is.finite(data[[y]])
  ggplot(
    data[finite, ],
    aes(x = .data[[x]], y = .data[[y]], color = .data[[color_variable]])
  ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(alpha = 0.65, size = 1.5) +
    scale_color_manual(values = color_values, drop = FALSE) +
    coord_equal() +
    theme_comparison() +
    labs(
      title = title,
      x = paste("EM", axis_label),
      y = paste("MixSQP", axis_label),
      color = if (color_variable == "annotation") "Annotation" else "Variant class"
    )
}

log_likelihood_plot <- comparison_plot(
  model_plot_data, "log_likelihood_em", "log_likelihood_mixsqp",
  "Achieved log likelihood", "log likelihood"
)
mutvar_plot <- comparison_plot(
  model_plot_data, "mutational_variance_em", "mutational_variance_mixsqp",
  "Mutational variance", "mutational variance"
)
penetrance_plot <- comparison_plot(
  model_plot_data, "effective_penetrance_em", "effective_penetrance_mixsqp",
  "Effective penetrance", "effective penetrance"
)
enrichment_plot <- comparison_plot(
  enrichment_plot_data, "enrichment_em", "enrichment_mixsqp",
  "Mutational-variance enrichment", "enrichment",
  color_variable = "annotation", color_values = annotation_colors
) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 1, 5, 20, 100)) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 5, 20, 100))

plots <- list(
  log_likelihood = log_likelihood_plot,
  mutational_variance = mutvar_plot,
  effective_penetrance = penetrance_plot,
  enrichment = enrichment_plot
)
for (plot_name in names(plots)) {
  ggsave(
    file.path(component_figure_dir,
              paste0("AutismOptimizerComparison_", plot_name, ".pdf")),
    plots[[plot_name]], width = 11, height = 4.5
  )
}

combined_plot <- wrap_plots(plots, ncol = 2, guides = "collect") +
  plot_annotation(
    title = "Autism model comparison: MixSQP versus EM",
    tag_levels = "A"
  ) &
  theme(legend.position = "bottom")
ggsave(
  file.path(output_figure_dir, "AutismOptimizerComparison.pdf"),
  combined_plot, width = 14, height = 10
)

cat(
  "Matched", nrow(model_comparison), "model/prevalence fits across",
  length(unique(model_comparison$model_name)), "model names.\n"
)
cat(
  "Matched", nrow(enrichment_comparison),
  "model/prevalence/annotation enrichment estimates.\n"
)
cat("Comparison tables and figures written under Final_Runs_July2026/outputs.\n")
