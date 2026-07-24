library(ggplot2)
library(data.table)
library(patchwork)
library(readr)

script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(script_file)) {
  script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
}
script_file <- normalizePath(script_file, mustWork = TRUE)
repo_dir <- dirname(dirname(script_file))
final_runs_dir <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT", unset = file.path(repo_dir, "analysis")), mustWork = FALSE)

input_dir <- file.path(final_runs_dir, "inputs")
constraint_dir <- file.path(input_dir, "constraint")
age_input_dir <- file.path(input_dir, "age_phenotypes")
pedigree_dir <- file.path(age_input_dir, "pedigrees")
variant_dir <- file.path(age_input_dir, "variants")
count_dir <- file.path(input_dir, "autism", "gene_counts")
mutation_rate_dir <- file.path(input_dir, "mutation_rates")

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

theme_bhr_gridlines <- function(){
  theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          axis.text=element_text(size=15, color = "black"),
          axis.title=element_text(size=15,  color = "black"),
          legend.text=element_text(size=15),
          legend.title = element_text(size=15),
          legend.position = "None",
          legend.direction = "horizontal",
          strip.text.x = element_text(size = 15),
          strip.background = element_rect(fill = "white"))

}


#Function from Kyle for fixing Gene_IDs
resolve_duplicate_names_prior_to_collecting_by_gene_id <- function(input_table) {
  input_table$Gene_ID[input_table$Gene == 'HIST1H4F'] = 'ENSG00000274618'
  input_table$LOEUF[input_table$Gene == 'HIST1H4F'] = NA
  input_table$LOEUF_bin[input_table$Gene == 'HIST1H4F'] = NA

  input_table$Gene_ID[input_table$Gene == 'ATXN7'] = 'ENSG00000163635'
  input_table$LOEUF[input_table$Gene == 'ATXN7'] = 0.328
  input_table$LOEUF_bin[input_table$Gene == 'ATXN7'] = 1

  input_table$Gene_ID[input_table$Gene == 'PRSS50'] = 'ENSG00000283706'
  input_table$LOEUF[input_table$Gene == 'PRSS50'] = 0.99
  input_table$LOEUF_bin[input_table$Gene == 'PRSS50'] = 5

  input_table$Gene_ID[input_table$Gene == 'CCDC39'] = 'ENSG00000284862'
  input_table$LOEUF[input_table$Gene == 'CCDC39'] = 0.754
  input_table$LOEUF_bin[input_table$Gene == 'CCDC39'] = 4

  input_table$Gene_ID[input_table$Gene == 'HSPA14'] = 'ENSG00000187522'
  input_table$LOEUF[input_table$Gene == 'HSPA14'] = 0.395
  input_table$LOEUF_bin[input_table$Gene == 'HSPA14'] = 1

  input_table$Gene_ID[input_table$Gene == 'IGF2'] = 'ENSG00000167244'
  input_table$LOEUF[input_table$Gene == 'IGF2'] = 1.131
  input_table$LOEUF_bin[input_table$Gene == 'IGF2'] = 6

  input_table$Gene_ID[input_table$Gene == 'PDE11A'] = 'ENSG00000128655'
  input_table$LOEUF[input_table$Gene == 'PDE11A'] = 1.404
  input_table$LOEUF_bin[input_table$Gene == 'PDE11A'] = 7

  # These occur across multiple chromosomes
  temp1 = subset(input_table, Gene %in% c('RF00017', 'RF00019'))
  temp2 = subset(input_table, !Gene %in% c('RF00017', 'RF00019'))

  if(nrow(temp1) > 0) {
    temp1$Gene = paste0(temp1$Gene, '_chr', temp1$Chrom)
    input_table = rbind(temp1, temp2)
  } else {
    input_table = temp2
  }

  # RF00019 requires further parsing; here take the one that VEP uses most often
  input_table$Gene_ID[input_table$Gene == 'RF00019_chr3'] = 'ENSG00000212392'
  input_table$LOEUF[input_table$Gene == 'RF00019_chr3'] = NA
  input_table$LOEUF_bin[input_table$Gene == 'RF00019_chr3'] = NA

  input_table$Gene_ID[input_table$Gene == 'RF00019_chr6'] = 'ENSG00000200314'
  input_table$LOEUF[input_table$Gene == 'RF00019_chr6'] = NA
  input_table$LOEUF_bin[input_table$Gene == 'RF00019_chr6'] = NA

  input_table$Gene_ID[input_table$Gene == 'RF00019_chr7'] = 'ENSG00000201913'
  input_table$LOEUF[input_table$Gene == 'RF00019_chr7'] = NA
  input_table$LOEUF_bin[input_table$Gene == 'RF00019_chr7'] = NA

  input_table$Gene_ID[input_table$Gene == 'RF00019_chr12'] = 'ENSG00000207176'
  input_table$LOEUF[input_table$Gene == 'RF00019_chr12'] = NA
  input_table$LOEUF_bin[input_table$Gene == 'RF00019_chr12'] = NA

  input_table$Gene_ID[input_table$Gene == 'RF00019_chr16'] = 'ENSG00000199668'
  input_table$LOEUF[input_table$Gene == 'RF00019_chr16'] = NA
  input_table$LOEUF_bin[input_table$Gene == 'RF00019_chr16'] = NA

  return(input_table)
}

gnomad_information <- read.table(
  file.path(constraint_dir, "gnomad.v4.1.constraint_metrics.tsv"),
  header = TRUE
)

age_phenotypes <- read.table(
  file.path(age_input_dir, "20260312_asc_age_phenotypes.tsv"),
  sep = "\t", header = TRUE
)
age_phenotypes_spark = age_phenotypes[age_phenotypes$cohort == 'spark' & !is.na(age_phenotypes$diagnosis_age),]
age_phenotypes_spark$diagnosis_age_Above6 = as.numeric(age_phenotypes_spark$diagnosis_age > 72)

ped_files = file.path(pedigree_dir, c("SPARK_iWES_v2_de_novo_fam_v1.1c_new_samples_only_2025-03-29_AgePhenoProband.txt",
              "SPARK_Pilot_GATK_published_fam_for_de_novo_calls_2025-03-29_AgePhenoProband.txt",
              "SPARK_WES1_GATK_published_fam_for_de_novo_calls_2025-03-29_AgePhenoProband.txt"))

variant_files = file.path(variant_dir, c("SPARK_iWES_v2_de_novo_calls_v1.1c_new_samples_only_2025-03-29_AgePhenoProband.txt",
                  "SPARK_Pilot_GATK_published_autosomal_and_updated_XY_de_novo_calls_2025-03-29_AgePhenoProband.txt",
                  "SPARK_WES1_GATK_published_autosomal_and_updated_XY_de_novo_calls_2025-03-29_AgePhenoProband.txt"))

gene_count_file <- file.path(
  count_dir,
  "SPARK_iWES_v2_de_novo_counts_by_gene_new_samples_only_new_mis_cats_2025-03-29.txt"
)

variant_data <- lapply(variant_files, read.table, header = TRUE, sep = "\t")
names(variant_data) <- basename(variant_files)



#get the PTV-PTV_NonIndel ratios
num_ptv_nonindels <- sapply(variant_data,
                            function(data) {
                              return(sum((data$isPTV | data$isOS) & !data$isIndel))
                            })

num_ptv <- sapply(variant_data,
                  function(data) {
                    return(sum((data$isPTV | data$isOS)))
                  })


scale_ptv_factor = sum(num_ptv)/sum(num_ptv_nonindels)

mutation_rate_table = read.table(file.path(mutation_rate_dir, "ASD_gene_table_w_bespoke_mutation_rates_2024-07-24.txt"),
                                 header = TRUE,
                                 sep = "\t")
library(SummarizedExperiment)

#get the number of proband de novos for each variant class

proband_ptv = sapply(variant_data,
                     function(data) {
                       return(sum((data$isPTV | data$isOS) & data$Role == "Proband"))
                     })

proband_ptv_indel = sapply(variant_data,
                           function(data) {
                             return(sum((data$isPTV | data$isOS) &data$isIndel & data$Role == "Proband"))
                           })

proband_mis2 = sapply(variant_data,
                      function(data) {
                        alphamissense_filter = as.numeric(!is.na(data$am_pathogenicity) & data$am_pathogenicity >= 0.97)
                        MPC2_filter = as.numeric(!is.na(data$MPC_v2) & data$MPC_v2 >= 2)

                        missense_filter = (alphamissense_filter + MPC2_filter) == 2
                        return(sum((data$isMis & !data$isOS & missense_filter) & data$Role == "Proband"))
                      })

proband_mis1 = sapply(variant_data,
                      function(data) {
                        alphamissense_filter = as.numeric(!is.na(data$am_pathogenicity) & data$am_pathogenicity >= 0.97)
                        MPC2_filter = as.numeric(!is.na(data$MPC_v2) & data$MPC_v2 >= 2)

                        missense_filter = (alphamissense_filter + MPC2_filter) == 1
                        return(sum((data$isMis& !data$isOS & missense_filter) & data$Role == "Proband"))
                      })

proband_mis0 = sapply(variant_data,
                      function(data) {
                        alphamissense_filter = as.numeric(!is.na(data$am_pathogenicity) & data$am_pathogenicity >= 0.97)
                        MPC2_filter = as.numeric(!is.na(data$MPC_v2) & data$MPC_v2 >= 2)

                        missense_filter = (alphamissense_filter + MPC2_filter) == 0
                        return(sum((data$isMis& !data$isOS & missense_filter) & data$Role == "Proband"))
                      })


proband_syn = sapply(variant_data,
                     function(data) {

                       return(sum((data$isSyn & !data$isOS) & data$Role == "Proband"))
                     })


sum(proband_ptv) + sum(proband_mis2) + sum(proband_mis1) + sum(proband_mis0)+sum(proband_syn)
#First, construct the ColData


studies = c("SPARK iWES v2",
            "SPARK Pilot",
            "SPARK WES1 GATK")


phase = c("New",
          "Fu2022",
          "Fu2022")

colData <- data.frame()

for (filenum in 1:length(ped_files)) {
  print(filenum)
  file = ped_files[filenum]
  dataset = read.table(file, header = TRUE, sep = "\t")

  summary_df = data.frame(phase = phase[filenum],
                          study = studies[filenum],
                          N_Proband = sum(dataset$Role == "Proband"),
                          N_Sib = sum(dataset$Role == "Sibling"),
                          N_diagnosis_age_Above6 = sum(dataset$diagnosis_age_Above6 == 1),
                          N_diagnosis_age_AtOrBelow6 = sum(dataset$diagnosis_age_Above6 == 0))

  colData = rbind(colData,summary_df)
}

rownames(colData) = colData$study

#Now, define the rowData
gene_count_data <- read.table(gene_count_file, header = TRUE, sep = "\t")
rowData = data.frame(gene_id = gene_count_data$Gene_ID,
                     gene = gene_count_data$Gene,
                     LOEUF = gene_count_data$LOEUF,
                     Chrom = gene_count_data$Chrom)
rowData = cbind(rowData, mutation_rate_table[,-c(1:5)])
rownames(rowData) = rowData$gene_id


#initialize the SummarizedExperiment

autism_counts = SummarizedExperiment(rowData = rowData,
                                     colData = colData)

#Now, get the counts

tally_function <- function(variantdata, gene_ids, variant_class, Role, Age) {
  variantdata_fixGeneID = resolve_duplicate_names_prior_to_collecting_by_gene_id(variantdata)
  if (variant_class == "PTV") {
    variantdata_fixGeneID_subset = variantdata_fixGeneID[(variantdata_fixGeneID$isPTV | variantdata_fixGeneID$isOS)  & !variantdata_fixGeneID$isIndel & variantdata_fixGeneID$Role == Role & variantdata_fixGeneID$diagnosis_age_Above6 == Age,]
    return(sapply(gene_ids,
                  function(id) {
                    sum(variantdata_fixGeneID_subset$Gene_ID == id)
                  }))

  } else if (grepl("Mis",variant_class)) {
    num_criteria = parse_number(variant_class)
    variantdata_fixGeneID_missense = variantdata_fixGeneID[variantdata_fixGeneID$isMis & !variantdata_fixGeneID$isOS & variantdata_fixGeneID$Role == Role & variantdata_fixGeneID$diagnosis_age_Above6 == Age, ]
    alphamissense_filter = as.numeric(!is.na(variantdata_fixGeneID_missense$am_pathogenicity) & variantdata_fixGeneID_missense$am_pathogenicity >= 0.97)
    MPC2_filter = as.numeric(!is.na(variantdata_fixGeneID_missense$MPC_v2) & variantdata_fixGeneID_missense$MPC_v2 >= 2)

    missense_filter = (alphamissense_filter + MPC2_filter) == num_criteria
    variantdata_fixGeneID_missense_filter = variantdata_fixGeneID_missense[missense_filter,]

    return(sapply(gene_ids,
                  function(id) {
                    sum(variantdata_fixGeneID_missense_filter$Gene_ID == id)
                  }))
  } else if (variant_class == "Syn") {
    variantdata_fixGeneID_subset = variantdata_fixGeneID[variantdata_fixGeneID$isSyn & !variantdata_fixGeneID$isOS & variantdata_fixGeneID$Role == Role & variantdata_fixGeneID$diagnosis_age_Above6 == Age,]
    return(sapply(gene_ids,
                  function(id) {
                    sum(variantdata_fixGeneID_subset$Gene_ID == id)
                  }))
  }
}




count_subset_names = c("PTV_Proband_Age1","PTV_Proband_Age0",
                       "Mis2_Proband_Age1","Mis2_Proband_Age0",
                       "Mis1_Proband_Age1","Mis1_Proband_Age0",
                       "Mis0_Proband_Age1","Mis0_Proband_Age0",
                       "Syn_Proband_Age1","Syn_Proband_Age0")

subset_variant_class = stringr::str_split_i(count_subset_names,"_",1)
subset_Role = stringr::str_split_i(count_subset_names,"_",2)
subset_Age = as.numeric(stringr::str_split_i(count_subset_names,"_Age",2))


for (subsetnum in 1:length(count_subset_names)) {
  print(count_subset_names[subsetnum])
  subset_count_matrix = matrix(data = NA, nrow = nrow(rowData), ncol = nrow(colData))
  rownames(subset_count_matrix) = rownames(rowData)
  colnames(subset_count_matrix) = rownames(colData)


  for (num_count in seq_along(variant_data)) {
    print(num_count)

    variantdata = variant_data[[num_count]]

    subset_count_matrix[,num_count] = tally_function(variantdata,rownames(subset_count_matrix),
                                                     subset_variant_class[subsetnum],
                                                     subset_Role[subsetnum],
                                                     subset_Age[subsetnum])
  }

  assays(autism_counts)[[count_subset_names[subsetnum]]] = subset_count_matrix
}

#subset to genes with complete mutation rate information
autism_counts = autism_counts[complete.cases(mutation_rate_table[,-c(1:5)]),]

#subset to genes with available LOEUF information
autism_counts = autism_counts[!is.na(rowData(autism_counts)$LOEUF),]

#subset to autosomal genes
autism_counts = autism_counts[!(rowData(autism_counts)$Chrom %in% c("X","Y")),]

#make a one-hot encoded LOEUF annotation matrix
quantiles <- quantile(rowData(autism_counts)$LOEUF, probs = seq(0, 1, by = 0.2))

# Use cut to map each value to its quantile range
quantile_assignments <- cut(rowData(autism_counts)$LOEUF, breaks = quantiles, include.lowest = TRUE,
                            labels = c(1:5))

for (bin in 1:5) {
  rowData(autism_counts)[paste0("LOEUF",bin)] =  as.numeric(quantile_assignments == bin)
}


#Get the posterior mutation rate correction factors based on gnomad counts
#First, let's get the mutation rate correction factors

gnomad_input_df <- data.frame(case_count = gnomad_information$syn.obs[gnomad_information$canonical == "true" & grepl("ENSG",gnomad_information$gene_id)],
                              expected_count = gnomad_information$syn.exp[gnomad_information$canonical == "true"& grepl("ENSG",gnomad_information$gene_id)],
                              gene_id = gnomad_information$gene_id[gnomad_information$canonical == "true"& grepl("ENSG",gnomad_information$gene_id)])
gnomad_input_df = gnomad_input_df[complete.cases(gnomad_input_df),]

mut_calibration_syn <- BurdenMLE_DN(gnomad_input_df,
                                     features = NULL,
                                     component_endpoints = seq(-2,2,length.out = 31),
                                     mutvar_est = FALSE,
                                     max_iter = 1000,
                                     null_sim = FALSE,
                                     bootstrap = FALSE,
                                     return_likelihood = TRUE,
                                     estimate_posteriors = TRUE,
                                     estimate_effective_penetrance = FALSE,
                                     optimizer = "mixsqp")

rowData(autism_counts)$PosteriorMuCorrectionFactor = mut_calibration_syn$posterior_gene_estimates$Posterior_Mean[match(rowData(autism_counts)$gene_id,
                                                                                                                       gnomad_input_df$gene_id)]
ggplot(mapping = aes(x = rowData(autism_counts)$PosteriorMuCorrectionFactor))+
  geom_histogram(color = "black", fill = "white")+
  theme_bhr_legend_gridlines()+
  labs(x = "Posterior Mean\nMutation Rate Correction", y = "Count")

autism_counts = autism_counts[!is.na(rowData(autism_counts)$PosteriorMuCorrectionFactor),]

#Define the "baseline model" features
#get the class-wise mutation rates
LOEUF = as.matrix(rowData(autism_counts)[,13:17])

mu_byclass = as.matrix(rowData(autism_counts)[,5:9])
cumulative_mu = rowSums(mu_byclass) * rowData(autism_counts)$PosteriorMuCorrectionFactor
cumulative_mu_bin <- cbind(as.numeric(cumulative_mu > quantile(cumulative_mu, 0.8)),
                           as.numeric(cumulative_mu <= quantile(cumulative_mu, 0.8)))


autism_features <- cbind(LOEUF[,1] * cumulative_mu_bin[,1],
                         LOEUF[,1] * cumulative_mu_bin[,2],
                         LOEUF[,2:5])
colnames(autism_features) = c("LOEUF1_mu1",
                              "LOEUF1_mu2",
                              "LOEUF2",
                              "LOEUF3",
                              "LOEUF4",
                              "LOEUF5")


autism_loop_vars = list()
# Only the two age strata used by the final visualization are retained.
autism_loop_vars$names = c(
  "SPARK Probands AgeAbove6 PTV", "SPARK Probands AgeAtOrBelow6 PTV",
  "SPARK Probands AgeAbove6 Mis2", "SPARK Probands AgeAtOrBelow6 Mis2",
  "SPARK Probands AgeAbove6 Mis1", "SPARK Probands AgeAtOrBelow6 Mis1",
  "SPARK Probands AgeAbove6 Mis0", "SPARK Probands AgeAtOrBelow6 Mis0",
  "SPARK Probands AgeAbove6 Syn", "SPARK Probands AgeAtOrBelow6 Syn"
)

autism_loop_vars$datasets = rep(
  list(autism_counts$phase %in% c("New", "Fu2022")), 10
)

autism_loop_vars$N_subset = rep(
  c(sum(autism_counts$N_diagnosis_age_Above6),
    sum(autism_counts$N_diagnosis_age_AtOrBelow6)),
  5
)

autism_loop_vars$subsets = list(
  "PTV_Proband_Age1", "PTV_Proband_Age0",
  "Mis2_Proband_Age1", "Mis2_Proband_Age0",
  "Mis1_Proband_Age1", "Mis1_Proband_Age0",
  "Mis0_Proband_Age1", "Mis0_Proband_Age0",
  "Syn_Proband_Age1", "Syn_Proband_Age0"
)

autism_loop_vars$mutation_rate = rep(
  c("mu_snp_PTV", "mu_snp_Mis2", "mu_snp_Mis1", "mu_snp_Mis0", "mu_snp_Syn"),
  each = 2
)

baseprev = 0.0276

autism_loop_vars$prevalences = rep(baseprev, 10)
autism_loop_vars$fit_model = rep(TRUE, 10)

loop_var_lengths <- vapply(autism_loop_vars, length, integer(1))
if (length(unique(loop_var_lengths)) != 1) {
  stop("Age-phenotype loop-variable fields have different lengths: ",
       paste(names(loop_var_lengths), loop_var_lengths, sep = "=", collapse = ", "))
}


# Retain the main-analysis prevalence reference while deriving the autism
# cohort sizes from the current pedigrees. DDD is not needed for age-specific
# counts, but remains part of the combined-cohort prevalence reference.
N_SPARK = sum(autism_counts$N_Proband[grepl("SPARK", autism_counts$study)])
N_ASC = sum(autism_counts$N_Proband[grepl("ASC", autism_counts$study)])
N_GeneDX = sum(autism_counts$N_Proband[grepl("GeneDx", autism_counts$study)])
N_DDD = 1194
N_total = N_SPARK + N_ASC + N_GeneDX + N_DDD

weight_SPARK = N_SPARK/N_total
weight_ASC = N_ASC/N_total
weight_GeneDX = N_GeneDX/N_total
weight_DDD = N_DDD/N_total

prev_weighted = (weight_SPARK * 0.0276) + (weight_ASC * 0.0276) + (weight_GeneDX * 0.01) + (weight_DDD * 0.01)

prev_factors = c(prev_weighted/0.0276, 0.01/0.0276, 1)


autism_data <- list(counts = autism_counts,
                    features = autism_features,
                    loop_vars = autism_loop_vars,
                    ptv_scale_factor = scale_ptv_factor,
                    baseprev = baseprev,
                    prev_factors = prev_factors)
