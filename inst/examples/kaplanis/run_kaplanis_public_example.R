library(BurdenMLEDN)
library(ggplot2)

kaplanis_url <- paste0(
  "https://media.springernature.com/original/springer-static/esm/",
  "art%3A10.1038%2Fs41586-020-2832-5/MediaObjects/",
  "41586_2020_2832_MOESM3_ESM.txt"
)

n_trios <- 31058
prevalence <- 0.01
n_boot <- 100

# In RStudio, set the working directory to the script's directory:
# Session > Set Working Directory > To Source File Location
output_dir <- getwd()
variant_file <- file.path(
  output_dir,
  "Kaplanis2020_SupplementaryTable1.txt"
)

download.file(kaplanis_url, variant_file, mode = "wb")

variants <- read.delim(variant_file, stringsAsFactors = FALSE)
ptv_consequences <- c(
  "frameshift_variant",
  "splice_acceptor_variant",
  "splice_donor_variant",
  "stop_gained"
)

variants$class <- ifelse(
  variants$consequence %in% ptv_consequences,
  "PTV",
  ifelse(variants$consequence == "synonymous_variant", "Synonymous", NA)
)
variants$is_indel <- nchar(variants$ref) != nchar(variants$alt)

# The bundled PTV rates cover SNVs, so scale them to include observed indels.
ptv <- variants$class == "PTV" & !is.na(variants$class)
ptv_scale <- sum(ptv) / sum(ptv & !variants$is_indel)

# The package reference supplies calibrated mutation rates and six LOEUF strata.
reference <- burdenmle_gene_reference()

gene_counts <- function(class_name) {
  counts <- table(variants$symbol[variants$class == class_name & !is.na(variants$class)])
  result <- as.integer(counts[match(reference$gene_symbol, names(counts))])
  result[is.na(result)] <- 0L
  result
}

features <- loeuf_features(reference$gene_id, reference)

fit_class <- function(class_name, mutation_rate, rate_scale = 1) {
  input <- data.frame(
    case_count = gene_counts(class_name),
    case_rate = mutation_rate * rate_scale,
    N = n_trios,
    row.names = reference$gene_id
  )

  fit <- BurdenMLE_DN(
    input,
    features = features,
    component_endpoints = seq(0, log(250), length.out = 10),
    prevalence = prevalence,
    n_boot = n_boot,
    bootstrap_seed = 2020,
    null_sim = FALSE,
    estimate_effective_penetrance = FALSE
  )

  ci <- as.numeric(fit$mutvar_output$mutvar_CI)
  data.frame(
    variant_class = class_name,
    observed_variants = sum(input$case_count),
    expected_variants = sum(2 * input$N * input$case_rate),
    mutational_variance = fit$mutvar_output$total_mutvar,
    lower_95_ci = ci[1],
    upper_95_ci = ci[2]
  )
}

results <- rbind(
  fit_class("PTV", reference$corrected_mu_snp_PTV, ptv_scale),
  fit_class("Synonymous", reference$corrected_mu_snp_Syn)
)

write.table(
  results,
  file.path(output_dir, "KaplanisPublicExample_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

figure <- ggplot(results, aes(variant_class, mutational_variance)) +
  geom_errorbar(
    aes(ymin = lower_95_ci, ymax = upper_95_ci),
    width = 0.12,
    linewidth = 0.7
  ) +
  geom_point(aes(fill = variant_class), shape = 21, size = 4, show.legend = FALSE) +
  scale_fill_manual(values = c(PTV = "#D95F5F", Synonymous = "#7A7A7A")) +
  scale_y_continuous(
    labels = function(x) sprintf("%.1f%%", 100 * x),
    limits = c(0, max(results$upper_95_ci) * 1.08),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = "Variant class", y = "Observed-scale mutational variance") +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA)
  )

ggsave(
  file.path(output_dir, "KaplanisPublicExample_mutational_variance.pdf"),
  figure,
  width = 6.4,
  height = 4.8
)

print(results)
