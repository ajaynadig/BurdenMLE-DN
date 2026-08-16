# Compare paternal age among SPARK probands with paternal age at U.S. births.
#
# SPARK father age is linked exactly as in age_sandbox.R: the three SPARK
# pedigree files are combined, phenotype data are matched by Sample, and the
# analysis is restricted to probands with a reported father_age_birth.

library(ggplot2)
library(patchwork)

script_file <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
)
script_file <- normalizePath(script_file, mustWork = TRUE)
repo_dir <- dirname(dirname(dirname(script_file)))
final_runs_dir <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT", unset = file.path(repo_dir, "analysis")), mustWork = FALSE)

pedigree_dir <- file.path(final_runs_dir, "inputs", "autism", "pedigrees")
phenotype_file <- file.path(
  final_runs_dir, "inputs", "age_phenotypes",
  "20260312_asc_age_phenotypes.tsv"
)
us_distribution_file <- file.path(
  final_runs_dir, "outputs", "tables", "supplementary", "paternal_age",
  "us_father_age_distribution_2024.csv"
)
table_dir <- file.path(final_runs_dir, "outputs", "tables", "supplementary", "paternal_age")
figure_dir <- file.path(final_runs_dir, "outputs", "figures", "supplementary")
diagnostic_figure_dir <- file.path(
  final_runs_dir, "outputs", "figures", "diagnostics", "paternal_age"
)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostic_figure_dir, recursive = TRUE, showWarnings = FALSE)

pedigree_files <- c(
  "SPARK_iWES_v2_de_novo_fam_v1.1c_new_samples_only_2025-03-29.txt",
  "SPARK_Pilot_GATK_published_fam_for_de_novo_calls_2025-03-29.txt",
  "SPARK_WES1_GATK_published_fam_for_de_novo_calls_2025-03-29.txt"
)

spark_pedigree <- do.call(
  rbind,
  lapply(file.path(pedigree_dir, pedigree_files), read.table, header = TRUE)
)
age_phenotypes <- read.delim(phenotype_file)
us_distribution <- read.csv(us_distribution_file)

if (anyDuplicated(spark_pedigree$Sample)) {
  stop("A Sample appears more than once across the combined SPARK pedigrees.")
}
if (anyDuplicated(age_phenotypes$Sample)) {
  stop("A Sample appears more than once in the age phenotype file.")
}

spark_probands <- spark_pedigree[spark_pedigree$Role == "Proband", ]
spark_probands$father_age_months <- age_phenotypes$father_age_birth[
  match(spark_probands$Sample, age_phenotypes$Sample)
]
number_all_probands <- nrow(spark_probands)
spark_probands <- spark_probands[is.finite(spark_probands$father_age_months), ]

# CDC/NCHS reports completed age in integer years. Floor the month-resolved
# SPARK ages so that summaries and plot bins use the same age definition.
spark_probands$father_age <- floor(spark_probands$father_age_months / 12)
if (any(spark_probands$father_age < 10 | spark_probands$father_age > 98)) {
  warning("SPARK contains father ages outside the CDC/NCHS reported range 10-98.")
}

weighted_quantile <- function(values, weights, probabilities) {
  ordering <- order(values)
  values <- values[ordering]
  weights <- weights[ordering]
  cumulative <- cumsum(weights) / sum(weights)
  vapply(probabilities, function(p) values[which(cumulative >= p)[1]], numeric(1))
}

summarize_ages <- function(ages, weights = rep(1, length(ages))) {
  n <- sum(weights)
  mean_age <- weighted.mean(ages, weights)
  data.frame(
    n = n,
    mean = mean_age,
    sd = sqrt(sum(weights * (ages - mean_age)^2) / (n - 1)),
    q1 = weighted_quantile(ages, weights, 0.25),
    median = weighted_quantile(ages, weights, 0.5),
    q3 = weighted_quantile(ages, weights, 0.75),
    percent_35_or_older = 100 * sum(weights[ages >= 35]) / n,
    percent_40_or_older = 100 * sum(weights[ages >= 40]) / n,
    percent_45_or_older = 100 * sum(weights[ages >= 45]) / n
  )
}

spark_summary <- summarize_ages(spark_probands$father_age)
us_summary <- summarize_ages(us_distribution$father_age, us_distribution$births)
summary_table <- rbind(
  cbind(group = "SPARK autistic probands", spark_summary),
  cbind(group = "U.S. births (2024)", us_summary)
)
summary_table$father_age_reported_percent <- c(
  100 * nrow(spark_probands) / number_all_probands,
  NA_real_
)

thresholds <- c(35, 40, 45)
threshold_plot_data <- do.call(rbind, lapply(thresholds, function(threshold) {
  x <- sum(spark_probands$father_age >= threshold)
  interval <- binom.test(x, nrow(spark_probands))$conf.int * 100
  us_percent <- 100 * sum(
    us_distribution$births[us_distribution$father_age >= threshold]
  ) / sum(us_distribution$births)
  rbind(
    data.frame(
      group = "SPARK autistic probands", threshold = threshold,
      percent = 100 * x / nrow(spark_probands), lower = interval[1],
      upper = interval[2]
    ),
    data.frame(
      group = "U.S. births (2024)", threshold = threshold,
      percent = us_percent, lower = us_percent, upper = us_percent
    )
  )
}))

mean_difference <- spark_summary$mean - us_summary$mean
mean_difference_se <- sqrt(
  spark_summary$sd^2 / spark_summary$n + us_summary$sd^2 / us_summary$n
)
welch_df <- mean_difference_se^4 / (
  (spark_summary$sd^2 / spark_summary$n)^2 / (spark_summary$n - 1) +
    (us_summary$sd^2 / us_summary$n)^2 / (us_summary$n - 1)
)
welch_t <- mean_difference / mean_difference_se
welch_p <- 2 * pt(-abs(welch_t), df = welch_df)
pooled_sd <- sqrt(
  ((spark_summary$n - 1) * spark_summary$sd^2 +
     (us_summary$n - 1) * us_summary$sd^2) /
    (spark_summary$n + us_summary$n - 2)
)
cohens_d <- mean_difference / pooled_sd

mean_test <- data.frame(
  test = "Welch two-sample t-test from sample summaries",
  estimate_spark = spark_summary$mean,
  estimate_us_2024 = us_summary$mean,
  mean_difference_years = mean_difference,
  ci_lower = mean_difference - qt(0.975, welch_df) * mean_difference_se,
  ci_upper = mean_difference + qt(0.975, welch_df) * mean_difference_se,
  t_statistic = welch_t,
  degrees_freedom = welch_df,
  p_value = welch_p,
  cohens_d = cohens_d
)

comparison_table <- rbind(
  data.frame(
    statistic = "Mean paternal age (years)",
    spark = spark_summary$mean, us_2024 = us_summary$mean,
    spark_minus_us = mean_difference,
    difference_ci_lower = mean_test$ci_lower,
    difference_ci_upper = mean_test$ci_upper
  ),
  do.call(rbind, lapply(thresholds, function(threshold) {
    rows <- threshold_plot_data[threshold_plot_data$threshold == threshold, ]
    spark_row <- rows[rows$group == "SPARK autistic probands", ]
    us_row <- rows[rows$group == "U.S. births (2024)", ]
    data.frame(
      statistic = paste0("Percent age ", threshold, " or older"),
      spark = spark_row$percent, us_2024 = us_row$percent,
      spark_minus_us = spark_row$percent - us_row$percent,
      difference_ci_lower = spark_row$lower - us_row$percent,
      difference_ci_upper = spark_row$upper - us_row$percent
    )
  }))
)

spark_age_counts <- as.data.frame(table(spark_probands$father_age))
names(spark_age_counts) <- c("father_age", "births")
spark_age_counts$father_age <- as.integer(as.character(spark_age_counts$father_age))
spark_age_counts$percent <- 100 * spark_age_counts$births / sum(spark_age_counts$births)
spark_age_counts$group <- "SPARK autistic probands"
us_plot_data <- transform(
  us_distribution,
  percent = percent_known_age,
  group = "U.S. births (2024)"
)
distribution_plot_data <- rbind(
  spark_age_counts[, c("father_age", "percent", "group")],
  us_plot_data[, c("father_age", "percent", "group")]
)

group_colors <- c(
  "SPARK autistic probands" = "#D95F02",
  "U.S. births (2024)" = "#377EB8"
)

distribution_plot <- ggplot(
  distribution_plot_data,
  aes(father_age, percent, color = group)
) +
  geom_line(linewidth = 1) +
  coord_cartesian(xlim = c(15, 65)) +
  scale_color_manual(values = group_colors) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.title.position = "plot"
  ) +
  labs(
    title = "A  Distribution of paternal age",
    subtitle = paste0(
      "SPARK: ", format(nrow(spark_probands), big.mark = ","),
      " probands with reported paternal age"
    ),
    x = "Paternal age at birth (completed years)",
    y = "Percent of births"
  )

threshold_plot_data$threshold_label <- factor(
  paste0("Age >=", threshold_plot_data$threshold),
  levels = paste0("Age >=", thresholds)
)
threshold_plot <- ggplot(
  threshold_plot_data,
  aes(threshold_label, percent, color = group, group = group)
) +
  geom_point(position = position_dodge(width = 0.35), size = 2.8) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.35), width = 0.12
  ) +
  scale_color_manual(values = group_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title.position = "plot"
  ) +
  labs(
    title = "B  Older paternal-age thresholds",
    subtitle = "SPARK error bars are exact 95% binomial confidence intervals",
    x = NULL,
    y = "Percent of births"
  )

comparison_figure <- distribution_plot / threshold_plot +
  plot_annotation(
    title = "Paternal age in SPARK and U.S. births",
    caption = paste0(
      "U.S. reference: 2024 CDC/NCHS natality records with reported father age. ",
      "SPARK ages were converted from months to completed years."
    )
  )

write.csv(
  summary_table,
  file.path(table_dir, "spark_vs_us_father_age_summary.csv"),
  row.names = FALSE
)
write.csv(
  threshold_plot_data[, c("group", "threshold", "percent", "lower", "upper")],
  file.path(table_dir, "spark_vs_us_father_age_thresholds.csv"),
  row.names = FALSE
)
write.csv(
  comparison_table,
  file.path(table_dir, "spark_vs_us_father_age_comparison.csv"),
  row.names = FALSE
)
write.csv(
  mean_test,
  file.path(table_dir, "spark_vs_us_father_age_mean_test.csv"),
  row.names = FALSE
)
ggsave(
  file.path(figure_dir, "SupplementaryFigure4_PaternalAge.pdf"),
  comparison_figure, width = 8, height = 8.5
)
ggsave(
  file.path(diagnostic_figure_dir, "SupplementaryFigure4_PaternalAge.png"),
  comparison_figure, width = 8, height = 8.5, dpi = 300
)

print(summary_table, row.names = FALSE)
print(mean_test, row.names = FALSE)
message(
  "SPARK father age was available for ", nrow(spark_probands), " of ",
  number_all_probands, " probands (",
  round(100 * nrow(spark_probands) / number_all_probands, 1), "%)."
)
message("Comparison tables and figure written under Final_Runs_July2026/outputs.")
