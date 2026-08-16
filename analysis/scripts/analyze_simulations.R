# Load required libraries and source files
library(ggplot2)
library(patchwork)
library(data.table)

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

normalize_optimizer <- function(x) {
  if (tolower(x) == "em") return("EM")
  if (tolower(x) == "mixsqp") return("mixsqp")
  stop("Optimizer must be EM or mixsqp.")
}

if ("--help" %in% args) {
  cat(
    "Options:\n",
    "  --optimizer <EM|mixsqp>          Results to analyze (default: mixsqp)\n",
    "  --compare-optimizer <EM|mixsqp>  Optional paired comparison\n",
    "  --run-date <label>                Output date label\n",
    "  --seed <integer>                  Simulation seed (default: 20260715)\n",
    sep = ""
  )
  quit(status = 0)
}

optimizer <- normalize_optimizer(get_arg("--optimizer", "mixsqp"))
compare_optimizer_arg <- get_arg("--compare-optimizer", "")
compare_optimizer <- if (compare_optimizer_arg == "") NULL else {
  normalize_optimizer(compare_optimizer_arg)
}
optimizer_label <- tolower(optimizer)
analysis_label <- optimizer_label
comparison_label <- if (is.null(compare_optimizer)) optimizer_label else {
  paste0(optimizer_label, "_vs_", tolower(compare_optimizer))
}

script_file <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
)
script_file <- normalizePath(script_file, mustWork = TRUE)
repo_dir <- dirname(dirname(dirname(script_file)))
final_runs_dir <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT", unset = file.path(repo_dir, "analysis")), mustWork = FALSE)
output_data_dir <- file.path(final_runs_dir, "outputs", "data")
main_figures_dir <- file.path(final_runs_dir, "outputs", "figures", "main")
supplementary_figures_dir <- file.path(final_runs_dir, "outputs", "figures", "supplementary")
diagnostic_figures_dir <- file.path(final_runs_dir, "outputs", "figures", "diagnostics", "simulations")
diagnostic_tables_dir <- file.path(final_runs_dir, "outputs", "tables", "diagnostics", "simulations")
supplementary_tables_dir <- file.path(final_runs_dir, "outputs", "tables", "supplementary")

for (directory in c(main_figures_dir, supplementary_figures_dir,
                    diagnostic_figures_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}
dir.create(diagnostic_tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supplementary_tables_dir, recursive = TRUE, showWarnings = FALSE)

run_date <- get_arg("--run-date", format(Sys.Date(), "%b%d_%y"))
if (!grepl("^[A-Za-z0-9_-]+$", run_date)) {
  stop("--run-date may contain only letters, numbers, underscores, and hyphens.")
}

seed_str <- get_arg("--seed", Sys.getenv("SIM_RNG_SEED"))
if (seed_str == "") seed_str <- "20260715"
rng_seed <- as.integer(seed_str)
if (is.na(rng_seed)) stop("--seed must be an integer.")

results_file_for_optimizer <- function(method) {
  file.path(
    output_data_dir,
    sprintf("simulation_results_%s_%s_seed%s.RData",
            tolower(method), run_date, rng_seed)
  )
}

simulation_results_file <- results_file_for_optimizer(optimizer)
effect_size_samples_file <- file.path(output_data_dir,
  sprintf("effect_size_samples_%s_seed%s.RData", run_date, rng_seed))

if(!file.exists(simulation_results_file)) stop(sprintf("Missing simulation results file: %s", simulation_results_file))
if(!file.exists(effect_size_samples_file)) stop(sprintf("Missing effect size samples file: %s", effect_size_samples_file))

# Load simulation results and effect size samples
load(simulation_results_file)
load(effect_size_samples_file)

theme_bhr_legend_gridlines <- function(){
  theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          axis.text=element_text(size=15, color = "black"),
          axis.title=element_text(size=15,  color = "black"),
          legend.text=element_text(size=15),
          legend.title = element_text(size=15),
          strip.text.x = element_text(size = 15),
          strip.background = element_rect(fill = "white"))
}

theme_bhr_legend <- function(){
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15,  color = "black"),
        legend.text=element_text(size=15),
        legend.title = element_text(size=15),
        strip.text.x = element_text(size = 15),
        strip.background = element_rect(fill = "white"))
}



# Function to create histogram plot from samples
create_effect_size_density_plot <- function(effect_size_df) {
  # Filter out null distribution
  effect_size_df <- subset(effect_size_df, distribution != "null")
  
  # Convert distribution to factor with specific order
  effect_size_df$distribution <- factor(effect_size_df$distribution,
                                      levels=c("half_uniform", "half_infinitesimal", "point_normal", "oligogenic"))
  
  ggplot(effect_size_df, aes(x = effect_size)) +
    geom_histogram(bins = 50, color = "black", fill = "white") +
    scale_y_log10() +  # Log scale y-axis
    theme_bw() +
    labs(x = "Log Rate Ratio", y = "Number of Genes") +
    facet_wrap(~distribution, scales = "free_x")  # Allow x scales to vary but keep y fixed
}

# Function to create summary dataframe from simulation results
create_summary_df <- function(results, optimizer_name = optimizer) {
  summary_df <- data.frame()
  
  for(dist in names(results)) {
    if(dist != "null") {  # Skip null simulations
      dist_results <- results[[dist]]
      
      # Extract metrics for each simulation
      get_metric <- function(x, name, default = NA_real_) {
        if (is.null(x[[name]])) default else x[[name]]
      }
      metrics <- data.frame(
        distribution = dist,
        simulation = seq_along(dist_results),
        N = sapply(dist_results, function(x) x$N),
        true_mutvar = sapply(dist_results, function(x) x$true_mutvar),
        est_mutvar = sapply(dist_results, function(x) x$estimated_mutvar),
        true_frac_cases = sapply(dist_results, function(x) x$true_frac_cases),
        est_frac_cases = sapply(dist_results, function(x) x$estimated_frac_cases),
        true_peneff = sapply(dist_results, function(x) x$true_peneff),
        est_peneff = sapply(dist_results, function(x) x$estimated_peneff),
        optimizer = sapply(dist_results, function(x) {
          get_metric(x, "optimizer", optimizer_name)
        }),
        log_likelihood = sapply(dist_results, function(x) {
          get_metric(x, "log_likelihood")
        }),
        log_likelihood_per_gene = sapply(dist_results, function(x) {
          get_metric(x, "log_likelihood_per_gene")
        }),
        total_fit_seconds = sapply(dist_results, function(x) {
          get_metric(x, "total_fit_seconds")
        }),
        optimizer_full_fit_seconds = sapply(dist_results, function(x) {
          get_metric(x, "optimizer_full_fit_seconds")
        }),
        optimizer_iterations = sapply(dist_results, function(x) {
          get_metric(x, "optimizer_iterations")
        }),
        optimizer_converged = sapply(dist_results, function(x) {
          get_metric(x, "optimizer_converged", NA)
        }),
        optimizer_status = sapply(dist_results, function(x) {
          get_metric(x, "optimizer_status", NA_character_)
        }),
        active_components_min = sapply(dist_results, function(x) {
          values <- get_metric(x, "active_components_by_stratum", NA_real_)
          min(values)
        }),
        active_components_mean = sapply(dist_results, function(x) {
          values <- get_metric(x, "active_components_by_stratum", NA_real_)
          mean(values)
        }),
        active_components_max = sapply(dist_results, function(x) {
          values <- get_metric(x, "active_components_by_stratum", NA_real_)
          max(values)
        })
      )
      metrics$simulation_within_N <- ave(
        metrics$simulation,
        metrics$N,
        FUN = seq_along
      )
      
      summary_df <- rbind(summary_df, metrics)
    }
  }
  
  # Convert distribution to factor with specific order
  summary_df$distribution <- factor(summary_df$distribution, 
                                  levels=c("half_uniform", "half_infinitesimal", "point_normal", "oligogenic"))
  
  # Convert N to factor for coloring
  summary_df$N <- factor(summary_df$N)
  
  return(summary_df)
}

# Create plots
create_simulation_plots <- function(summary_df,
                                    effect_size_df,
                                    output_label = analysis_label) {
  # Create distribution mapping for facet labels
  dist_labels <- c(
    "half_uniform" = "Half-Uniform",
    "half_infinitesimal" = "Half-Infinitesimal",
    "point_normal" = "Point-Normal",
    "oligogenic" = "Oligogenic"
  )

  # Create histogram plot
  histogram_plot <- create_effect_size_density_plot(effect_size_df) +
    facet_wrap(~distribution, scales="free", labeller = labeller(distribution = dist_labels))
  
  # Helper function to set axis limits based on data
  set_equal_limits <- function(data, true_col, est_col) {
    min_val <- min(c(data[[true_col]], data[[est_col]]))
    max_val <- max(c(data[[true_col]], data[[est_col]]))
    # Add small buffer
    range <- max_val - min_val
    buffer <- 0.05 * range
    return(c(min_val - buffer, max_val + buffer))
  }
  
  # Create limit dataframes for each distribution and N
  mutvar_limits_df <- do.call(rbind, lapply(unique(summary_df$distribution), function(dist) {
    do.call(rbind, lapply(unique(summary_df$N), function(n) {
      dist_data <- subset(summary_df, distribution == dist & N == n)
      limits <- set_equal_limits(dist_data, "true_mutvar", "est_mutvar")
      data.frame(
        distribution = dist,
        N = n,
        true_mutvar = limits,
        est_mutvar = limits
      )
    }))
  }))
  
  frac_limits_df <- do.call(rbind, lapply(unique(summary_df$distribution), function(dist) {
    do.call(rbind, lapply(unique(summary_df$N), function(n) {
      dist_data <- subset(summary_df, distribution == dist & N == n)
      limits <- set_equal_limits(dist_data, "true_frac_cases", "est_frac_cases")
      data.frame(
        distribution = dist,
        N = n,
        true_frac_cases = limits,
        est_frac_cases = limits
      )
    }))
  }))
  
  peneff_limits_df <- do.call(rbind, lapply(unique(summary_df$distribution), function(dist) {
    do.call(rbind, lapply(unique(summary_df$N), function(n) {
      dist_data <- subset(summary_df, distribution == dist & N == n)
      limits <- set_equal_limits(dist_data, "true_peneff", "est_peneff")
      data.frame(
        distribution = dist,
        N = n,
        true_peneff = limits,
        est_peneff = limits
      )
    }))
  }))
  
  # Define color palette and distribution labels
  sample_colors <- c(
    "1000" = "#3B528BFF",    # Deep blue
    "5000" = "#21908CFF",    # Teal
    "15000" = "#5DC863FF",   # Green
    "38680" = "#FDE725FF",   # Yellow
    "50000" = "#FEAF77FF"    # Light orange - softer endpoint that follows naturally
  )


  
  # 1. Mutational variance estimation plot
  mutvar_plot <- ggplot(summary_df, aes(x=true_mutvar, y=est_mutvar, color=N)) +
    geom_point(alpha=0.5) +
    geom_smooth(method="lm", se=FALSE) +
    geom_abline(intercept=0, slope=1, linetype="dashed", color="black") +
    geom_blank(data=mutvar_limits_df) +
    scale_color_manual(values = sample_colors, guide = "none") +
    theme_bw() +
    labs(x="True Mutational Variance", y="Estimated Mutational Variance") +
    facet_wrap(~distribution, scales="free", labeller = labeller(distribution = dist_labels)) +
    scale_x_continuous(expand = c(0.02, 0)) +
    scale_y_continuous(expand = c(0.02, 0))
  
  # 2. Fraction of cases plot
  frac_cases_plot <- ggplot(summary_df, aes(x=true_frac_cases, y=est_frac_cases, color=N)) +
    geom_point(alpha=0.5) +
    geom_smooth(method="lm", se=FALSE) +
    geom_abline(intercept=0, slope=1, linetype="dashed", color="black") +
    geom_blank(data=frac_limits_df) +
    scale_color_manual(values = sample_colors, name = "Number of Trios") +
    theme_bw() +
    labs(x="True Fraction of Cases (RR>5)", y="Estimated Fraction of Cases (RR>5)") +
    facet_wrap(~distribution, scales="free", labeller = labeller(distribution = dist_labels)) +
    scale_x_continuous(expand = c(0.02, 0)) +
    scale_y_continuous(expand = c(0.02, 0)) +
    theme(legend.position = "bottom")
  
  # 3. Effective penetrance plot
  peneff_plot <- ggplot(summary_df, aes(x=true_peneff, y=est_peneff, color=N)) +
    geom_point(alpha=0.5) +
    geom_smooth(method="lm", se=FALSE) +
    geom_abline(intercept=0, slope=1, linetype="dashed", color="black") +
    geom_blank(data=peneff_limits_df) +
    scale_color_manual(values = sample_colors, guide = "none") +
    theme_bw() +
    labs(x="True Effective Penetrance", y="Estimated Effective Penetrance") +
    facet_wrap(~distribution, scales="free", labeller = labeller(distribution = dist_labels)) +
    scale_x_continuous(expand = c(0.02, 0)) +
    scale_y_continuous(expand = c(0.02, 0))
  
  # Combine plots with histogram plot on top and add lettering
  combined_plot <- (histogram_plot + labs(tag = "A")) / 
    ((mutvar_plot + labs(tag = "B")) + 
     (frac_cases_plot + labs(tag = "C")) + 
     (peneff_plot + labs(tag = "D"))) +
    plot_layout(heights = c(1, 1)) &
    theme(plot.tag = element_text(size = 20, face = "bold"))
  
  return(combined_plot)
}

# Create Figure 1 with histograms and mutational variance estimates
create_figure_1 <- function(effect_size_df, summary_df) {
  # Define color palette
  sample_colors <- c(
    "1000" = "#3B528BFF",    # Deep blue
    "5000" = "#21908CFF",    # Teal
    "15000" = "#5DC863FF",   # Green
    "38680" = "#FDE725FF",   # Yellow
    "50000" = "#FEAF77FF"    # Light orange - softer endpoint that follows naturally
  )

  # Create distribution mapping
  dist_labels <- c(
    "half_infinitesimal" = "Half-Infinitesimal",
    "point_normal" = "Point-Normal",
    "oligogenic" = "Oligogenic"
  )
  
  # Filter data for specific distributions
  distributions_to_keep <- names(dist_labels)
  
  # Filter effect size data
  filtered_effect_size_df <- subset(effect_size_df, 
    distribution %in% distributions_to_keep)
  filtered_effect_size_df$distribution <- factor(filtered_effect_size_df$distribution,
    levels = distributions_to_keep)
  
  # Filter summary data
  filtered_summary_df <- subset(summary_df, 
    distribution %in% distributions_to_keep)
  filtered_summary_df$distribution <- factor(filtered_summary_df$distribution,
    levels = distributions_to_keep)

  # Per-panel annotation describing which parameters vary across simulations
  histogram_annotations <- data.frame(
    distribution = factor(distributions_to_keep, levels = distributions_to_keep),
    label_body = c(
      "Effect Size Spread",
      "Proportion Non-Null\nMean Effect Size\nEffect Size Spread",
      "Proportion Non-Null\nEffect Size"
    )
  )
  
  # Facet-specific annotation anchors (manually tuned for final figure layout)
  annotation_layout <- data.frame(
    distribution = factor(distributions_to_keep, levels = distributions_to_keep),
    x_pos = c(2.4, 2.4, 2.4),
    y_header = c(9500, 9500, 9500),
    y_body = c(5200, 5200, 5200)
  )
  histogram_annotations <- merge(histogram_annotations, annotation_layout, by = "distribution", all.x = TRUE)
  histogram_annotations$distribution <- factor(histogram_annotations$distribution, levels = distributions_to_keep)
  
  # Create histograms
  histogram_row <- ggplot(filtered_effect_size_df, aes(x = effect_size)) +
    geom_histogram(binwidth = 0.2, color = "black", fill = "white") +
    geom_text(
      data = histogram_annotations,
      aes(x = x_pos, y = y_header, label = "Parameters Varied"),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 3.2,
      lineheight = 1.05,
      fontface = "bold"
    ) +
    geom_text(
      data = histogram_annotations,
      aes(x = x_pos, y = y_body, label = label_body),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.8,
      lineheight = 1.12
    ) +
    scale_y_log10() +
    theme_bhr_legend() +
    labs(x = "Log Rate Ratio", y = "Number of Genes") +
    facet_wrap(~distribution, nrow = 1, labeller = labeller(distribution = dist_labels)) +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(5.5, 12, 5.5, 5.5))
  
  # Create mutational variance estimates plot
  mut_var_row <- ggplot(filtered_summary_df, aes(x = true_mutvar, y = est_mutvar, color = N)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    scale_color_manual(values = sample_colors) +
    theme_bhr_legend_gridlines() +
    labs(x = "True Mutational Variance", y = "Estimated\nMutational Variance", color = "Number of Trios") +
    facet_wrap(~distribution, nrow = 1, labeller = labeller(distribution = function(x) rep("", length(x)))) +
    theme(legend.position = "bottom",
          axis.line = element_line(color = "black", linewidth = 1),
          strip.background = element_blank(),
          strip.text = element_blank())
  
  # Combine plots
  combined_plot <- (histogram_row / mut_var_row) +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(tag_levels = 'A')&
  theme(plot.tag = element_text(size = 20, face = "bold"))
  
  return(combined_plot)
}

# Create weight comparison plot for half-uniform simulations
create_weight_comparison_plot <- function(results,
                                          output_label = analysis_label) {
  # Extract true and estimated weights from half-uniform simulations
  half_uniform_results <- results$half_uniform
  
  weight_df <- data.frame()
  
  for(i in seq_along(half_uniform_results)) {
    sim <- half_uniform_results[[i]]
    
    # Get estimated weights from model
    est_weights <- sim$model_weights  # Changed from sim$model$delta[1,]
    
    # Create dataframe for this simulation
    sim_df <- data.frame(
      N = sim$N,
      simulation = i,
      component = 1:length(est_weights),
      estimated_weight = est_weights,
      true_weight = sim$params$weights
    )
    
    weight_df <- rbind(weight_df, sim_df)
  }
  
  # Convert N to factor for coloring
  weight_df$N <- factor(weight_df$N)
  
  # Define a more balanced color palette
  # Using viridis-inspired colors with smoother progression
  sample_colors <- c(
    "1000" = "#3B528BFF",    # Deep blue
    "5000" = "#21908CFF",    # Teal
    "15000" = "#5DC863FF",   # Green
    "38680" = "#FDE725FF",   # Yellow
    "50000" = "#FEAF77FF"    # Light orange - softer endpoint that follows naturally
  )
  
  # Create scatterplot of true vs estimated weights
  weight_plot <- ggplot(weight_df, aes(x=true_weight, y=estimated_weight, color=N)) +
    geom_point(alpha=0.5) +
    geom_smooth(method="lm", se=FALSE) +
    geom_abline(intercept=0, slope=1, linetype="dashed", color="black") +
    scale_color_manual(values = sample_colors) +
    theme_bw() +
    labs(x="True Weight", y="Estimated Weight", 
         title="True vs Estimated Weights for Half-Uniform Distribution") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "right")
  
  return(weight_plot)
}

create_optimizer_comparison_plot <- function(paired_summary,
                                             primary_method,
                                             comparison_method) {
  primary_label <- tolower(primary_method)
  comparison_label <- tolower(comparison_method)
  sample_sizes <- sort(unique(as.numeric(as.character(paired_summary$N))))
  distribution_colors <- c(
    "half_uniform" = "#3B528BFF",
    "half_infinitesimal" = "#21908CFF",
    "point_normal" = "#5DC863FF",
    "oligogenic" = "#FEAF77FF"
  )

  runtime_df <- rbind(
    data.frame(
      N = paired_summary$N,
      distribution = paired_summary$distribution,
      optimizer = primary_method,
      seconds = paired_summary[[paste0("total_fit_seconds_", primary_label)]]
    ),
    data.frame(
      N = paired_summary$N,
      distribution = paired_summary$distribution,
      optimizer = comparison_method,
      seconds = paired_summary[[paste0("total_fit_seconds_", comparison_label)]]
    )
  )
  runtime_df$optimizer <- factor(
    runtime_df$optimizer,
    levels = c(comparison_method, primary_method)
  )

  plot_list <- list()
  for (row_index in seq_along(sample_sizes)) {
    sample_size <- sample_sizes[row_index]
    paired_N <- paired_summary[as.numeric(as.character(paired_summary$N)) == sample_size, ]
    runtime_N <- runtime_df[as.numeric(as.character(runtime_df$N)) == sample_size, ]

    runtime_plot <- ggplot(runtime_N, aes(x = optimizer, y = seconds, fill = optimizer)) +
      geom_boxplot(outlier.alpha = 0.15) +
      guides(fill = "none") +
      theme_bhr_legend_gridlines() +
      labs(
        title = if (row_index == 1) "Runtime" else NULL,
        subtitle = paste0("N = ", format(sample_size, big.mark = ",")),
        x = NULL,
        y = "Total fit time (seconds)"
      )

    likelihood_plot <- ggplot(
      paired_N,
      aes(
        x = .data[[paste0("log_likelihood_", comparison_label)]],
        y = .data[[paste0("log_likelihood_", primary_label)]],
        color = distribution
      )
    ) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      geom_point(alpha = 0.45, size = 1) +
      scale_color_manual(values = distribution_colors) +
      coord_equal() +
      theme_bhr_legend_gridlines() +
      labs(
        title = if (row_index == 1) "Achieved log likelihood" else NULL,
        x = paste(comparison_method, "log likelihood"),
        y = paste(primary_method, "log likelihood"),
        color = "Simulation distribution"
      )

    mutvar_plot <- ggplot(
      paired_N,
      aes(
        x = .data[[paste0("est_mutvar_", comparison_label)]],
        y = .data[[paste0("est_mutvar_", primary_label)]],
        color = distribution
      )
    ) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      geom_point(alpha = 0.45, size = 1) +
      scale_color_manual(values = distribution_colors) +
      coord_equal() +
      theme_bhr_legend_gridlines() +
      labs(
        title = if (row_index == 1) "Estimated mutational variance" else NULL,
        x = paste(comparison_method, "estimate"),
        y = paste(primary_method, "estimate"),
        color = "Simulation distribution"
      )

    plot_list <- c(plot_list, list(runtime_plot, likelihood_plot, mutvar_plot))
  }

  wrap_plots(plot_list, ncol = 3, guides = "collect") +
    plot_annotation(
      title = paste(primary_method, "and", comparison_method,
                    "simulation-fit comparison")
    ) &
    theme(legend.position = "bottom")
}

# Run analysis
summary_df <- create_summary_df(results, optimizer)
SupplementaryTable1 <- summary_df[, c(
  "distribution", "N", "true_mutvar", "est_mutvar",
  "true_frac_cases", "est_frac_cases", "true_peneff", "est_peneff"
)]
names(SupplementaryTable1) <- c(
  "Distribution",
  "N",
  "True Mutational Variance",
  "Estimated Mutational Variances",
  "True Fraction of Cases with RR>5 Variant",
  "Estimated Fraction of Cases with RR>5 Variant",
  "True Effective Penetrance",
  "Estimated Effective Penetrance"
)
write.table(SupplementaryTable1,
            file.path(diagnostic_tables_dir,
                      sprintf("SimulationEstimates_%s.tsv", toupper(optimizer_label))),
            sep = "\t", quote = FALSE, row.names = FALSE)
if (optimizer_label == "mixsqp") {
  write.table(
    SupplementaryTable1,
    file.path(
      supplementary_tables_dir,
      "SupplementaryTable1_SimulationResults.tsv"
    ),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

write.table(
  summary_df,
  file.path(
    diagnostic_tables_dir,
    sprintf("SimulationDiagnostics_%s.tsv", toupper(optimizer_label))
  ),
  sep = "\t", quote = FALSE, row.names = FALSE
)

if (!is.null(compare_optimizer)) {
  comparison_file <- results_file_for_optimizer(compare_optimizer)
  if (!file.exists(comparison_file)) {
    stop("Missing comparison results file: ", comparison_file)
  }
  comparison_environment <- new.env(parent = emptyenv())
  load(comparison_file, envir = comparison_environment)
  comparison_summary <- create_summary_df(
    comparison_environment$results,
    compare_optimizer
  )

  pair_columns <- c("distribution", "N", "simulation_within_N")
  truth_columns <- c("true_mutvar", "true_frac_cases", "true_peneff")
  paired_summary <- merge(
    summary_df,
    comparison_summary,
    by = pair_columns,
    suffixes = c(paste0("_", optimizer_label),
                 paste0("_", tolower(compare_optimizer))),
    sort = FALSE
  )
  if (nrow(paired_summary) != nrow(summary_df) ||
      nrow(paired_summary) != nrow(comparison_summary)) {
    stop("Optimizer result files do not contain the same simulation pairs.")
  }
  for (truth_column in truth_columns) {
    first <- paired_summary[[paste0(truth_column, "_", optimizer_label)]]
    second <- paired_summary[[paste0(truth_column, "_", tolower(compare_optimizer))]]
    if (!isTRUE(all.equal(first, second, tolerance = 0))) {
      stop(
        "Paired simulations do not have identical ", truth_column,
        ". Confirm the seed and simulation order."
      )
    }
  }

  primary_suffix <- paste0("_", optimizer_label)
  comparison_suffix <- paste0("_", tolower(compare_optimizer))
  paired_summary$estimated_mutvar_difference <-
    paired_summary[[paste0("est_mutvar", primary_suffix)]] -
    paired_summary[[paste0("est_mutvar", comparison_suffix)]]
  paired_summary$log_likelihood_difference <-
    paired_summary[[paste0("log_likelihood", primary_suffix)]] -
    paired_summary[[paste0("log_likelihood", comparison_suffix)]]
  paired_summary$total_runtime_speedup <-
    paired_summary[[paste0("total_fit_seconds", comparison_suffix)]] /
    paired_summary[[paste0("total_fit_seconds", primary_suffix)]]
  paired_summary$optimizer_runtime_speedup <-
    paired_summary[[paste0("optimizer_full_fit_seconds", comparison_suffix)]] /
    paired_summary[[paste0("optimizer_full_fit_seconds", primary_suffix)]]

  write.table(
    paired_summary,
    file.path(
      diagnostic_tables_dir,
      "SimulationOptimizerComparison.tsv"
    ),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  optimizer_comparison_plot <- create_optimizer_comparison_plot(
    paired_summary,
    optimizer,
    compare_optimizer
  )
  ggsave(
    file.path(
      supplementary_figures_dir,
      "SupplementaryFigure16_OptimizerComparison.pdf"
    ),
    optimizer_comparison_plot,
    width = 15,
    height = 20
  )
}

simulation_plots <- create_simulation_plots(
  summary_df,
  effect_size_df,
  optimizer_label
)
weight_plot <- create_weight_comparison_plot(results, optimizer_label)

# Print summary statistics
print("Summary Statistics by Distribution:")
for(dist in unique(summary_df$distribution)) {
  dist_data <- subset(summary_df, distribution == dist)
  
  cat("\nDistribution:", dist, "\n")
  cat("Mean MutVar Bias:", mean(dist_data$est_mutvar - dist_data$true_mutvar), "\n")
  cat("MutVar Bias SD:", sd(dist_data$est_mutvar - dist_data$true_mutvar), "\n")
  cat("Mean Frac Cases Bias:", mean(dist_data$est_frac_cases - dist_data$true_frac_cases), "\n")
  cat("Frac Cases Bias SD:", sd(dist_data$est_frac_cases - dist_data$true_frac_cases), "\n")
}

# Create and save all plots
summary_df <- create_summary_df(results, optimizer)
simulation_plots <- create_simulation_plots(
  summary_df,
  effect_size_df,
  optimizer_label
)
weight_plot <- create_weight_comparison_plot(results, optimizer_label)
figure_1 <- create_figure_1(effect_size_df, summary_df)

# Save plots with date in filenames
#Saving full results
print("Saving simulation results...")
full_results_file <- if (tolower(analysis_label) == "mixsqp") {
  file.path(supplementary_figures_dir, "SupplementaryFigure1_FullSimulationResults.pdf")
} else {
  file.path(diagnostic_figures_dir, sprintf("FullSimulationResults_%s.pdf", toupper(analysis_label)))
}
ggsave(full_results_file,
       simulation_plots, width=20 * (4/5), height=15* (2/3))
# Saving weight comparison plot
print("Saving weight comparison plot...")
ggsave(file.path(diagnostic_figures_dir,
                 sprintf("WeightComparison_%s.pdf", toupper(analysis_label))),
       weight_plot, width=8, height=6)
# Saving figure 1
print("Saving figure 1...")
figure1_file <- if (tolower(analysis_label) == "mixsqp") {
  file.path(main_figures_dir, "Figure1.pdf")
} else {
  file.path(diagnostic_figures_dir, sprintf("Figure1_%s.pdf", toupper(analysis_label)))
}
ggsave(figure1_file,
       figure_1, width=12, height=8)

if (!is.null(compare_optimizer)) {
  comparison_optimizer_label <- tolower(compare_optimizer)
  comparison_simulation_plots <- create_simulation_plots(
    comparison_summary,
    effect_size_df,
    comparison_optimizer_label
  )
  comparison_weight_plot <- create_weight_comparison_plot(
    comparison_environment$results,
    comparison_optimizer_label
  )
  comparison_figure_1 <- create_figure_1(
    effect_size_df,
    comparison_summary
  )

  ggsave(
    file.path(
      diagnostic_figures_dir,
      sprintf("FullSimulationResults_%s.pdf", toupper(comparison_optimizer_label))
    ),
    comparison_simulation_plots,
    width = 20 * (4 / 5), height = 15 * (2 / 3)
  )
  ggsave(
    file.path(
      diagnostic_figures_dir,
      sprintf("WeightComparison_%s.pdf", toupper(comparison_optimizer_label))
    ),
    comparison_weight_plot,
    width = 8, height = 6
  )
  ggsave(
    file.path(
      diagnostic_figures_dir,
      sprintf("Figure1_%s.pdf", toupper(comparison_optimizer_label))
    ),
    comparison_figure_1,
    width = 12, height = 8
  )
}
