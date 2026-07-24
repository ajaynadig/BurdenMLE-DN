# Describe paternal age at U.S. births using the 2024 CDC/NCHS Natality
# Public Use File. The compressed download is about 227 MB, but this script
# streams father age directly from the archive and does not extract the full
# multi-gigabyte fixed-width file.

library(ggplot2)

year <- 2024
source_url <- paste0(
  "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/DVS/natality/",
  "Nat", year, "us.zip"
)
documentation_url <- paste0(
  "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/",
  "Dataset_Documentation/DVS/natality/UserGuide", year, ".pdf"
)

script_file <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
)
script_file <- normalizePath(script_file, mustWork = TRUE)
repo_dir <- dirname(dirname(script_file))
final_runs_dir <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT", unset = file.path(repo_dir, "analysis")), mustWork = FALSE)

input_dir <- file.path(final_runs_dir, "inputs", "demographics")
table_dir <- file.path(final_runs_dir, "outputs", "tables", "supplementary", "paternal_age")
figure_dir <- file.path(final_runs_dir, "outputs", "figures", "diagnostics", "paternal_age")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

zip_file <- file.path(input_dir, paste0("Nat", year, "us.zip"))
if (!file.exists(zip_file)) {
  message("Downloading CDC/NCHS natality file (approximately 227 MB)...")
  download.file(source_url, zip_file, mode = "wb", quiet = FALSE)
}

archive_contents <- unzip(zip_file, list = TRUE)
data_member <- archive_contents$Name[
  which.max(archive_contents$Length)
]
message("Reading father age from archive member: ", data_member)

# FAGECOMB occupies fixed-width positions 147-148 in the 2024 user guide.
# Values 10-98 are reported ages; 99 and ages under 10 are not stated.
age_counts <- integer(100)
names(age_counts) <- as.character(0:99)
birth_count <- 0L
connection <- pipe(
  paste("unzip -p", shQuote(zip_file), shQuote(data_member)),
  open = "r"
)
on.exit(close(connection), add = TRUE)

repeat {
  lines <- readLines(connection, n = 100000L, warn = FALSE)
  if (length(lines) == 0) break
  # RESTATUS at position 104: codes 1-3 are residents of a U.S. state or DC;
  # code 4 identifies foreign residents whose births occurred in the U.S.
  residence_status <- suppressWarnings(as.integer(substr(lines, 104, 104)))
  lines <- lines[residence_status %in% 1:3]
  father_age <- suppressWarnings(as.integer(substr(lines, 147, 148)))
  father_age[is.na(father_age)] <- 99L
  age_counts <- age_counts + tabulate(father_age + 1L, nbins = 100L)
  birth_count <- birth_count + length(lines)
  if (birth_count %% 1000000L == 0L) {
    message(format(birth_count, big.mark = ","), " birth records processed")
  }
}
close(connection)
on.exit(NULL, add = FALSE)

known_ages <- 10:98
known_counts <- age_counts[as.character(known_ages)]
known_births <- sum(known_counts)
unknown_births <- birth_count - known_births

weighted_quantile <- function(values, weights, probabilities) {
  cumulative <- cumsum(weights) / sum(weights)
  vapply(probabilities, function(p) values[which(cumulative >= p)[1]], numeric(1))
}

mean_age <- weighted.mean(known_ages, known_counts)
sd_age <- sqrt(
  sum(known_counts * (known_ages - mean_age)^2) / (known_births - 1)
)
quantiles <- weighted_quantile(
  known_ages, known_counts, c(0.025, 0.25, 0.5, 0.75, 0.975)
)

summary_table <- data.frame(
  statistic = c(
    "Total U.S. births", "Births with reported father age",
    "Births without reported father age", "Percent father age missing",
    "Mean father age", "SD father age", "2.5th percentile", "25th percentile",
    "Median father age", "75th percentile", "97.5th percentile",
    "Percent age 35 or older", "Percent age 40 or older",
    "Percent age 45 or older"
  ),
  value = c(
    birth_count, known_births, unknown_births, 100 * unknown_births / birth_count,
    mean_age, sd_age, quantiles,
    100 * sum(known_counts[known_ages >= 35]) / known_births,
    100 * sum(known_counts[known_ages >= 40]) / known_births,
    100 * sum(known_counts[known_ages >= 45]) / known_births
  )
)

distribution <- data.frame(
  father_age = known_ages,
  births = as.integer(known_counts),
  percent_known_age = 100 * known_counts / known_births,
  percent_all_births = 100 * known_counts / birth_count
)

write.csv(
  summary_table,
  file.path(table_dir, paste0("us_father_age_summary_", year, ".csv")),
  row.names = FALSE
)
write.csv(
  distribution,
  file.path(table_dir, paste0("us_father_age_distribution_", year, ".csv")),
  row.names = FALSE
)
writeLines(
  c(
    paste("Data source:", source_url),
    paste("Documentation:", documentation_url),
    "Father age variable: FAGECOMB, fixed-width positions 147-148.",
    paste0(
      "Important limitation: father age is missing disproportionately for ",
      "births to unmarried mothers; statistics exclude records with unknown age."
    )
  ),
  file.path(table_dir, paste0("us_father_age_source_notes_", year, ".txt"))
)

age_distribution_plot <- ggplot(distribution, aes(father_age, percent_known_age)) +
  geom_col(width = 0.9, fill = "#2C7FB8") +
  geom_vline(xintercept = mean_age, linetype = "dashed", linewidth = 0.8) +
  annotate(
    "text", x = mean_age + 1, y = Inf,
    label = paste0("Mean = ", round(mean_age, 1)),
    hjust = 0, vjust = 1.5
  ) +
  coord_cartesian(xlim = c(15, 65)) +
  theme_bw(base_size = 12) +
  labs(
    title = paste("Father age at U.S. births,", year),
    subtitle = "Births with reported father age",
    x = "Father age (years)", y = "Percent of births"
  )

age_cdf <- transform(
  distribution,
  cumulative_percent = 100 * cumsum(births) / sum(births)
)
age_cdf_plot <- ggplot(age_cdf, aes(father_age, cumulative_percent)) +
  geom_step(linewidth = 0.9, color = "#D95F0E") +
  geom_hline(yintercept = c(25, 50, 75), linetype = "dotted") +
  coord_cartesian(xlim = c(15, 65)) +
  theme_bw(base_size = 12) +
  labs(
    title = paste("Cumulative distribution of father age,", year),
    subtitle = "Births with reported father age",
    x = "Father age (years)", y = "Cumulative percent"
  )

ggsave(
  file.path(figure_dir, paste0("us_father_age_distribution_", year, ".pdf")),
  age_distribution_plot, width = 8, height = 5
)
ggsave(
  file.path(figure_dir, paste0("us_father_age_cdf_", year, ".pdf")),
  age_cdf_plot, width = 8, height = 5
)

print(summary_table, row.names = FALSE)
message("Tables and figures written under Final_Runs_July2026/outputs.")


# SPARK paternal age and mutation-rate sensitivity analysis ------------------
#
# Quantify how quickly de novo mutation counts increase with paternal age in
# SPARK itself, then use that slope to assess the effect of the SPARK-versus-
# U.S. paternal-age difference on mutational variance. Synonymous variants are
# used because their counts provide a comparatively neutral mutation-rate
# proxy. The regression adjusts for SPARK sequencing dataset and child sex.

study_input_dir <- file.path(final_runs_dir, "inputs")
pedigree_dir <- file.path(study_input_dir, "autism", "pedigrees")
variant_dir <- file.path(study_input_dir, "autism", "variants")
age_phenotype_file <- file.path(
  study_input_dir, "age_phenotypes", "20260312_asc_age_phenotypes.tsv"
)

spark_pedigree_files <- file.path(
  pedigree_dir,
  c(
    "SPARK_iWES_v2_de_novo_fam_v1.1c_new_samples_only_2025-03-29.txt",
    "SPARK_Pilot_GATK_published_fam_for_de_novo_calls_2025-03-29.txt",
    "SPARK_WES1_GATK_published_fam_for_de_novo_calls_2025-03-29.txt"
  )
)
spark_variant_files <- file.path(
  variant_dir,
  c(
    "SPARK_iWES_v2_de_novo_calls_v1.1c_new_samples_only_2025-03-29.txt",
    "SPARK_Pilot_GATK_published_autosomal_and_updated_XY_de_novo_calls_2025-03-29.txt",
    "SPARK_WES1_GATK_published_autosomal_and_updated_XY_de_novo_calls_2025-03-29.txt"
  )
)

spark_pedigree <- data.table::rbindlist(lapply(spark_pedigree_files, data.table::fread))
age_phenotypes <- data.table::fread(age_phenotype_file)
spark_probands <- spark_pedigree[Role == "Proband"]
spark_probands[, father_age_months := age_phenotypes$father_age_birth[
  match(Sample, age_phenotypes$Sample)
]]
spark_probands <- spark_probands[is.finite(father_age_months)]
spark_probands[, father_age_years := father_age_months / 12]
spark_probands[, father_age_completed := floor(father_age_years)]

spark_variants <- data.table::rbindlist(lapply(
  spark_variant_files,
  function(path) {
    data.table::fread(path, select = c("Sample", "Role", "Simplified_csq"))
  }
))
synonymous_counts <- spark_variants[
  Role == "Proband" & Simplified_csq == "synonymous_variant",
  .(synonymous_count = .N),
  by = Sample
]
spark_probands[, synonymous_count := synonymous_counts$synonymous_count[
  match(Sample, synonymous_counts$Sample)
]]
spark_probands[is.na(synonymous_count), synonymous_count := 0L]

# A log-link makes the paternal-age coefficient directly interpretable as the
# proportional change in mutation rate per year. Dataset adjustment controls
# for differences among the iWES, Pilot, and WES1 call sets.
mutation_age_model <- glm(
  synonymous_count ~ father_age_years + Dataset + Sex,
  family = poisson(link = "log"),
  data = spark_probands
)
age_coefficient <- coef(mutation_age_model)["father_age_years"]
age_coefficient_se <- sqrt(vcov(mutation_age_model)[
  "father_age_years", "father_age_years"
])

# Use completed years for both datasets because FAGECOMB in the CDC/NCHS file
# reports completed paternal age. This prevents a roughly half-year artifact
# from comparing month-resolved SPARK ages with integer U.S. ages.
spark_mean_completed_age <- mean(spark_probands$father_age_completed)
us_mean_completed_age <- weighted.mean(distribution$father_age, distribution$births)
mean_age_difference <- spark_mean_completed_age - us_mean_completed_age
mutation_rate_ratio <- exp(age_coefficient * mean_age_difference)
mutation_rate_ratio_ci <- exp(
  (age_coefficient + c(-1, 1) * 1.96 * age_coefficient_se) *
    mean_age_difference
)

mutation_sensitivity <- data.frame(
  statistic = c(
    "SPARK probands with paternal age",
    "Mean SPARK paternal age (completed years)",
    "Mean U.S. paternal age (completed years)",
    "SPARK minus U.S. mean age (years)",
    "Proportional mutation-rate increase per paternal year",
    "Expected SPARK/U.S. mutation-rate ratio",
    "Expected proportional mutational-variance bias"
  ),
  estimate = c(
    nrow(spark_probands), spark_mean_completed_age, us_mean_completed_age,
    mean_age_difference, exp(age_coefficient) - 1,
    mutation_rate_ratio, mutation_rate_ratio - 1
  ),
  lower_95 = c(
    NA, NA, NA, NA,
    exp(age_coefficient - 1.96 * age_coefficient_se) - 1,
    mutation_rate_ratio_ci[1], mutation_rate_ratio_ci[1] - 1
  ),
  upper_95 = c(
    NA, NA, NA, NA,
    exp(age_coefficient + 1.96 * age_coefficient_se) - 1,
    mutation_rate_ratio_ci[2], mutation_rate_ratio_ci[2] - 1
  )
)
write.csv(
  mutation_sensitivity,
  file.path(table_dir, "spark_paternal_age_mutation_rate_sensitivity.csv"),
  row.names = FALSE
)

# Apply the proportional correction to the main MixSQP PTV and Mis2 estimates
# when that model output is present. Model indices are matched by name because
# the streamlined model list may change order.
model_files <- list.files(
  file.path(final_runs_dir, "outputs", "data"),
  pattern = "^models_autism_mixsqp_.*\\.Rdata$",
  full.names = TRUE
)
if (length(model_files) > 0) {
  model_file <- model_files[which.max(file.info(model_files)$mtime)]
  model_environment <- new.env(parent = emptyenv())
  load(model_file, envir = model_environment)

  model_names <- model_environment$autism_data$loop_vars$names
  target_names <- c("Combined Probands PTV", "Combined Probands Mis2")
  target_indices <- match(target_names, model_names)
  if (!anyNA(target_indices)) {
    reported_mutvar <- vapply(
      target_indices,
      function(index) {
        model_environment$BurdenMLE_DN_models_autism[[1]][[index]]$mutvar_output$total_mutvar
      },
      numeric(1)
    )
    mutvar_sensitivity <- data.frame(
      variant_class = c(target_names, "Combined PTV + Mis2"),
      reported_mutvar = c(reported_mutvar, sum(reported_mutvar)),
      paternal_age_adjusted_mutvar = c(
        reported_mutvar / mutation_rate_ratio,
        sum(reported_mutvar) / mutation_rate_ratio
      )
    )
    mutvar_sensitivity$absolute_change <-
      mutvar_sensitivity$paternal_age_adjusted_mutvar -
      mutvar_sensitivity$reported_mutvar
    mutvar_sensitivity$percent_change <- 100 * (
      mutvar_sensitivity$paternal_age_adjusted_mutvar /
        mutvar_sensitivity$reported_mutvar - 1
    )
    write.csv(
      mutvar_sensitivity,
      file.path(table_dir, "spark_paternal_age_mutvar_sensitivity.csv"),
      row.names = FALSE
    )
  }
}

print(mutation_sensitivity, row.names = FALSE)
message(
  "The ", round(mean_age_difference, 3),
  "-year paternal-age difference implies an approximately ",
  round(100 * (mutation_rate_ratio - 1), 2),
  "% proportional mutation-rate and mutational-variance difference."
)
