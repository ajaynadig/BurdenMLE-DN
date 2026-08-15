library(ggplot2)
library(data.table)
library(patchwork)
library(readr)

# Locate paths from this script so the workflow does not depend on setwd().
script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(script_file)) {
  script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
}
script_file <- normalizePath(script_file, mustWork = TRUE)
repo_dir <- dirname(dirname(script_file))
final_runs_dir <- normalizePath(Sys.getenv("BURDENMLEDN_ANALYSIS_ROOT", unset = file.path(repo_dir, "analysis")), mustWork = FALSE)

input_dir <- file.path(final_runs_dir, "inputs")
constraint_dir <- file.path(input_dir, "constraint")
autism_input_dir <- file.path(input_dir, "autism")
pedigree_dir <- file.path(autism_input_dir, "pedigrees")
variant_dir <- file.path(autism_input_dir, "variants")
count_dir <- file.path(autism_input_dir, "gene_counts")
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

ped_files = file.path(pedigree_dir, c("SPARK_iWES_v2_de_novo_fam_v1.1c_new_samples_only_2025-03-29.txt",
              "ASC_B17_B21_de_novos_v2.2_ped_2025-03-29.txt",
              "GeneDx_ASC_de_novos_GRCh38_v3.4_ped_HPO_terms_removed_2025-07-05.txt",
              "ASC_v17_published_fam_for_de_novo_calls_2025-03-29.txt",
              "ASC_B15_B16_published_fam_for_de_novo_calls_2025-03-29.txt",
              "SPARK_Pilot_GATK_published_fam_for_de_novo_calls_2025-03-29.txt",
              "SPARK_WES1_GATK_published_fam_for_de_novo_calls_2025-03-29.txt",
              "Kaplanis_DDD_ASD_ped_1194_probands_2025-04-09.txt"))

# Retained as an input for provenance and for defining the common gene rows.
gene_count_file <- file.path(
  count_dir,
  "SPARK_iWES_v2_de_novo_counts_by_gene_new_samples_only_new_mis_cats_2025-03-29.txt"
)

variant_files = file.path(variant_dir, c("SPARK_iWES_v2_de_novo_calls_v1.1c_new_samples_only_2025-03-29.txt",
                  "ASC_B17_B21_de_novos_v2.2_calls_2025-03-29.txt",
                  "GeneDx_ASC_de_novos_GRCh38_v3.4_calls_2025-07-05.txt",
                  "ASC_v17_published_autosomal_and_updated_XY_de_novo_calls_2025-03-29.txt",
                  "ASC_B15_B16_published_autosomal_and_updated_XY_de_novo_calls_2025-03-29.txt",
                  "SPARK_Pilot_GATK_published_autosomal_and_updated_XY_de_novo_calls_2025-03-29.txt",
                  "SPARK_WES1_GATK_published_autosomal_and_updated_XY_de_novo_calls_2025-03-29.txt",
                  "Kaplanis_DDD_ASD_de_novos_GRCh38_2025-04-09.txt"))

# Each variant table is reused for several summaries and all assay strata.
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
            "ASC B17-B21",
            "GeneDx (GRCh38)",
            "ASC v17",
            "ASC B15-B16",
            "SPARK Pilot",
            "SPARK WES1 GATK",
            "DDD ASD")


phase = c("New",
          "New",
          "New",
          "Fu2022",
          "Fu2022",
          "Fu2022",
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
                          N_MaleProband = sum(dataset$Sex == "Male" & dataset$Role == "Proband"),
                          N_FemaleProband = sum(dataset$Sex == "Female"& dataset$Role == "Proband"),
                          N_DDID = sum(dataset$DDID[dataset$Role == "Proband"]),
                          N_MaleDDID = sum(dataset$Sex == "Male" & dataset$Role == "Proband" & dataset$DDID == 1),
                          N_MaleNonDDID = sum(dataset$Sex == "Male" & dataset$Role == "Proband" & dataset$DDID == 0),
                          N_FemaleDDID = sum(dataset$Sex == "Female" & dataset$Role == "Proband" & dataset$DDID == 1),
                          N_FemaleNonDDID = sum(dataset$Sex == "Female" & dataset$Role == "Proband" & dataset$DDID == 0))

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

tally_function <- function(variantdata, gene_ids, variant_class, Role, Sex, DDID_status) {
  variantdata_fixGeneID = resolve_duplicate_names_prior_to_collecting_by_gene_id(variantdata)
  if (variant_class == "PTV") {
    variantdata_fixGeneID_subset = variantdata_fixGeneID[(variantdata_fixGeneID$isPTV | variantdata_fixGeneID$isOS)  & !variantdata_fixGeneID$isIndel & variantdata_fixGeneID$Role == Role & variantdata_fixGeneID$Sex == Sex & variantdata_fixGeneID$DDID == DDID_status,]
    return(sapply(gene_ids,
                  function(id) {
                    sum(variantdata_fixGeneID_subset$Gene_ID == id)
                  }))

  } else if (grepl("Mis",variant_class)) {
    num_criteria = parse_number(variant_class)
    variantdata_fixGeneID_missense = variantdata_fixGeneID[variantdata_fixGeneID$isMis & !variantdata_fixGeneID$isOS & variantdata_fixGeneID$Role == Role & variantdata_fixGeneID$Sex == Sex & variantdata_fixGeneID$DDID == DDID_status, ]
    alphamissense_filter = as.numeric(!is.na(variantdata_fixGeneID_missense$am_pathogenicity) & variantdata_fixGeneID_missense$am_pathogenicity >= 0.97)
    MPC2_filter = as.numeric(!is.na(variantdata_fixGeneID_missense$MPC_v2) & variantdata_fixGeneID_missense$MPC_v2 >= 2)

    missense_filter = (alphamissense_filter + MPC2_filter) == num_criteria
    variantdata_fixGeneID_missense_filter = variantdata_fixGeneID_missense[missense_filter,]

    return(sapply(gene_ids,
                  function(id) {
                    sum(variantdata_fixGeneID_missense_filter$Gene_ID == id)
                  }))
  } else if (variant_class == "Syn") {
    variantdata_fixGeneID_subset = variantdata_fixGeneID[variantdata_fixGeneID$isSyn & !variantdata_fixGeneID$isOS & variantdata_fixGeneID$Role == Role & variantdata_fixGeneID$Sex == Sex & variantdata_fixGeneID$DDID == DDID_status,]
    return(sapply(gene_ids,
                  function(id) {
                    sum(variantdata_fixGeneID_subset$Gene_ID == id)
                  }))
  }
}




count_subset_names = c("PTV_Proband_Male_DDID1","PTV_Proband_Female_DDID1","PTV_Proband_Male_DDID0","PTV_Proband_Female_DDID0",
                       "PTV_Sibling_Male_DDID1","PTV_Sibling_Female_DDID1","PTV_Sibling_Male_DDID0","PTV_Sibling_Female_DDID0",
                       "Mis2_Proband_Male_DDID1","Mis2_Proband_Female_DDID1","Mis2_Proband_Male_DDID0","Mis2_Proband_Female_DDID0",
                       "Mis2_Sibling_Male_DDID1","Mis2_Sibling_Female_DDID1","Mis2_Sibling_Male_DDID0","Mis2_Sibling_Female_DDID0",
                       "Mis1_Proband_Male_DDID1","Mis1_Proband_Female_DDID1","Mis1_Proband_Male_DDID0","Mis1_Proband_Female_DDID0",
                       "Mis1_Sibling_Male_DDID1","Mis1_Sibling_Female_DDID1","Mis1_Sibling_Male_DDID0","Mis1_Sibling_Female_DDID0",
                       "Mis0_Proband_Male_DDID1","Mis0_Proband_Female_DDID1","Mis0_Proband_Male_DDID0","Mis0_Proband_Female_DDID0",
                       "Mis0_Sibling_Male_DDID1","Mis0_Sibling_Female_DDID1","Mis0_Sibling_Male_DDID0","Mis0_Sibling_Female_DDID0",
                       "Syn_Proband_Male_DDID1","Syn_Proband_Female_DDID1","Syn_Proband_Male_DDID0","Syn_Proband_Female_DDID0",
                       "Syn_Sibling_Male_DDID1","Syn_Sibling_Female_DDID1","Syn_Sibling_Male_DDID0","Syn_Sibling_Female_DDID0")

subset_variant_class = stringr::str_split_i(count_subset_names,"_",1)
subset_Role = stringr::str_split_i(count_subset_names,"_",2)
subset_Sex = stringr::str_split_i(count_subset_names,"_",3)
subset_DDID = readr::parse_number(stringr::str_split_i(count_subset_names,"_",4))


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
                                                     subset_Sex[subsetnum],
                                                     subset_DDID[subsetnum])
  }

  assays(autism_counts)[[count_subset_names[subsetnum]]] = subset_count_matrix
}

#Let's do a sanity check with the spark iWES2
sparkiWES2_counts = gene_count_data

sanitycheck_subsets = names(sparkiWES2_counts)[-c(1:5)]


for (subsetnum in 1:length(sanitycheck_subsets)) {
  print(sanitycheck_subsets[subsetnum])
  tallied_counts = (assays(autism_counts)[[paste0(sanitycheck_subsets[subsetnum],"_DDID1")]] + assays(autism_counts)[[paste0(sanitycheck_subsets[subsetnum],"_DDID0")]])[,1]
  print(table(tallied_counts == sparkiWES2_counts[,sanitycheck_subsets[subsetnum]]))
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
rownames(gnomad_input_df) <- gnomad_input_df$gene_id

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


#autism_counts = autism_counts[!is.na(rowData(autism_counts)$PosteriorMuCorrectionFactor_v2),]
autism_loop_vars = list()
# Define looping variables for the main BurdenMLE-DN loop.
autism_loop_vars$names = c("New Probands PTV",
                           "New Siblings PTV",
                           "Fu2022 Probands PTV",
                           "Fu2022 Siblings PTV",
                           "Combined Probands PTV",
                           "Combined Siblings PTV",
                           "New Probands Mis2",
                           "New Siblings Mis2",
                           "Fu2022 Probands Mis2",
                           "Fu2022 Siblings Mis2",
                           "Combined Probands Mis2",
                           "Combined Siblings Mis2",
                           "New Probands Mis1",
                           "New Siblings Mis1",
                           "Fu2022 Probands Mis1",
                           "Fu2022 Siblings Mis1",
                           "Combined Probands Mis1",
                           "Combined Siblings Mis1",
                           "New Probands Mis0",
                           "New Siblings Mis0",
                           "Fu2022 Probands Mis0",
                           "Fu2022 Siblings Mis0",
                           "Combined Probands Mis0",
                           "Combined Siblings Mis0",
                           "New Probands Syn",
                           "New Siblings Syn",
                           "Fu2022 Probands Syn",
                           "Fu2022 Siblings Syn",
                           "Combined Probands Syn",
                           "Combined Siblings Syn",
                           "Combined ASC PTV",
                           "Male ASC PTV",
                           "Female ASC PTV",
                           "Combined SPARK PTV",
                           "Male SPARK PTV",
                           "Female SPARK PTV",
                           "Combined GeneDx PTV",
                           "Male GeneDx PTV",
                           "Female GeneDx PTV",
                           "Combined Male PTV",
                           "Combined Female PTV",
                           "Combined Male Mis2",
                           "Combined Female Mis2",
                           "Combined DDID PTV",
                           "Combined Non-DDID PTV",
                           "Combined DDID Mis2",
                           "Combined Non-DDID Mis2",
                           "Combined DDID Mis1",
                           "Combined Non-DDID Mis1",
                           "Combined DDID Mis0",
                           "Combined Non-DDID Mis0",
                           "Combined DDID Syn",
                           "Combined Non-DDID Syn",
                           "SPARK + ASC PTV",
                           "SPARK + ASC Mis2",
                           "SPARK + ASC Mis1",
                           "SPARK + ASC Mis0",
                           "SPARK + ASC Syn",
                           "Male DDID ASC PTV",
                           "Male DDID ASC Mis2",
                           "Female DDID ASC PTV",
                           "Female DDID ASC Mis2",
                           "Male Non-DDID ASC PTV",
                           "Male Non-DDID ASC Mis2",
                           "Female Non-DDID ASC PTV",
                           "Female Non-DDID ASC Mis2",
                           "Male DDID SPARK PTV",
                           "Male DDID SPARK Mis2",
                           "Female DDID SPARK PTV",
                           "Female DDID SPARK Mis2",
                           "Male Non-DDID SPARK PTV",
                           "Male Non-DDID SPARK Mis2",
                           "Female Non-DDID SPARK PTV",
                           "Female Non-DDID SPARK Mis2",
                           "Male DDID GeneDx PTV",
                           "Male DDID GeneDx Mis2",
                           "Female DDID GeneDx PTV",
                           "Female DDID GeneDx Mis2",
                           "Male Non-DDID GeneDx PTV",
                           "Male Non-DDID GeneDx Mis2",
                           "Female Non-DDID GeneDx PTV",
                           "Female Non-DDID GeneDx Mis2",
                           "Male DDID PTV",
                           "Male Non-DDID PTV",
                           "Female DDID PTV",
                           "Female Non-DDID PTV",
                           "Male DDID Mis2",
                           "Male Non-DDID Mis2",
                           "Female DDID Mis2",
                           "Female Non-DDID Mis2")




autism_loop_vars$datasets = c(rep(list((autism_counts$phase == "New"),
                                       (autism_counts$phase == "New"),
                                       (autism_counts$phase == "Fu2022"),
                                       (autism_counts$phase == "Fu2022"),
                                       (autism_counts$phase %in% c("New","Fu2022")),
                                       (autism_counts$phase %in% c("New","Fu2022"))),
                                  5),
                              list(grepl("ASC",autism_counts$study),
                                   grepl("ASC",autism_counts$study),
                                   grepl("ASC",autism_counts$study),
                                   grepl("SPARK",autism_counts$study),
                                   grepl("SPARK",autism_counts$study),
                                   grepl("SPARK",autism_counts$study),
                                   grepl("GeneDx", autism_counts$study),
                                   grepl("GeneDx", autism_counts$study),
                                   grepl("GeneDx", autism_counts$study),
                                   (autism_counts$phase %in% c("New","Fu2022")),
                                   (autism_counts$phase %in% c("New","Fu2022")),
                                   (autism_counts$phase %in% c("New","Fu2022")),
                                   (autism_counts$phase %in% c("New","Fu2022")),
                                   (autism_counts$phase %in% c("New","Fu2022")),
                                   (autism_counts$phase %in% c("New","Fu2022")),
                                   (autism_counts$phase %in% c("New","Fu2022")),
                                   (autism_counts$phase %in% c("New","Fu2022")),
                                   (autism_counts$phase %in% c("New","Fu2022")),
                                   (autism_counts$phase %in% c("New","Fu2022")),
                                   (autism_counts$phase %in% c("New","Fu2022")),
                                   (autism_counts$phase %in% c("New","Fu2022")),
                                   (autism_counts$phase %in% c("New","Fu2022")),
                                   (autism_counts$phase %in% c("New","Fu2022")),
                                   (grepl("ASC",autism_counts$study) | grepl("SPARK",autism_counts$study)),
                                   (grepl("ASC",autism_counts$study) | grepl("SPARK",autism_counts$study)),
                                   (grepl("ASC",autism_counts$study) | grepl("SPARK",autism_counts$study)),
                                   (grepl("ASC",autism_counts$study) | grepl("SPARK",autism_counts$study)),
                                   (grepl("ASC",autism_counts$study) | grepl("SPARK",autism_counts$study)),
                                   grepl("ASC",autism_counts$study),
                                   grepl("ASC",autism_counts$study),
                                   grepl("ASC",autism_counts$study),
                                   grepl("ASC",autism_counts$study),
                                   grepl("ASC",autism_counts$study),
                                   grepl("ASC",autism_counts$study),
                                   grepl("ASC",autism_counts$study),
                                   grepl("ASC",autism_counts$study),
                                   grepl("SPARK",autism_counts$study),
                                   grepl("SPARK",autism_counts$study),
                                   grepl("SPARK",autism_counts$study),
                                   grepl("SPARK",autism_counts$study),
                                   grepl("SPARK",autism_counts$study),
                                   grepl("SPARK",autism_counts$study),
                                   grepl("SPARK",autism_counts$study),
                                   grepl("SPARK",autism_counts$study),
                                   grepl("GeneDx",autism_counts$study),
                                   grepl("GeneDx",autism_counts$study),
                                   grepl("GeneDx",autism_counts$study),
                                   grepl("GeneDx",autism_counts$study),
                                   grepl("GeneDx",autism_counts$study),
                                   grepl("GeneDx",autism_counts$study),
                                   grepl("GeneDx",autism_counts$study),
                                   grepl("GeneDx",autism_counts$study),
                                   (autism_counts$phase %in% c("New","Fu2022")),
                                   (autism_counts$phase %in% c("New","Fu2022")),
                                   (autism_counts$phase %in% c("New","Fu2022")),
                                   (autism_counts$phase %in% c("New","Fu2022")),
                                   (autism_counts$phase %in% c("New","Fu2022")),
                                   (autism_counts$phase %in% c("New","Fu2022")),
                                   (autism_counts$phase %in% c("New","Fu2022")),
                                   (autism_counts$phase %in% c("New","Fu2022"))))




autism_loop_vars$N_subset = c(rep(c(sum(autism_counts$N_Proband[autism_counts$phase == "New"]),
                                    sum(autism_counts$N_Sib[autism_counts$phase == "New"]),
                                    sum(autism_counts$N_Proband[autism_counts$phase == "Fu2022"]),
                                    sum(autism_counts$N_Sib[autism_counts$phase == "Fu2022"]),
                                    sum(autism_counts$N_Proband),
                                    sum(autism_counts$N_Sib)),
                                  5),
                              sum(autism_counts$N_Proband[grepl("ASC",autism_counts$study)]),
                              sum(autism_counts$N_MaleProband[grepl("ASC",autism_counts$study)]),
                              sum(autism_counts$N_FemaleProband[grepl("ASC",autism_counts$study)]),
                              sum(autism_counts$N_Proband[grepl("SPARK",autism_counts$study)]),
                              sum(autism_counts$N_MaleProband[grepl("SPARK",autism_counts$study)]),
                              sum(autism_counts$N_FemaleProband[grepl("SPARK",autism_counts$study)]),
                              sum(autism_counts$N_Proband[grepl("GeneDx",autism_counts$study)]),
                              sum(autism_counts$N_MaleProband[grepl("GeneDx",autism_counts$study)]),
                              sum(autism_counts$N_FemaleProband[grepl("GeneDx",autism_counts$study)]),
                              sum(autism_counts$N_MaleProband),
                              sum(autism_counts$N_FemaleProband),
                              sum(autism_counts$N_MaleProband),
                              sum(autism_counts$N_FemaleProband),
                              sum(autism_counts$N_DDID),
                              sum(autism_counts$N_Proband - autism_counts$N_DDID),
                              sum(autism_counts$N_DDID),
                              sum(autism_counts$N_Proband - autism_counts$N_DDID),
                              sum(autism_counts$N_DDID),
                              sum(autism_counts$N_Proband - autism_counts$N_DDID),
                              sum(autism_counts$N_DDID),
                              sum(autism_counts$N_Proband - autism_counts$N_DDID),
                              sum(autism_counts$N_DDID),
                              sum(autism_counts$N_Proband - autism_counts$N_DDID),
                              sum(autism_counts$N_Proband[(grepl("ASC",autism_counts$study) | grepl("SPARK",autism_counts$study))]),
                              sum(autism_counts$N_Proband[(grepl("ASC",autism_counts$study) | grepl("SPARK",autism_counts$study))]),
                              sum(autism_counts$N_Proband[(grepl("ASC",autism_counts$study) | grepl("SPARK",autism_counts$study))]),
                              sum(autism_counts$N_Proband[(grepl("ASC",autism_counts$study) | grepl("SPARK",autism_counts$study))]),
                              sum(autism_counts$N_Proband[(grepl("ASC",autism_counts$study) | grepl("SPARK",autism_counts$study))]),
                              sum(autism_counts$N_MaleDDID[grepl("ASC",autism_counts$study)]),
                              sum(autism_counts$N_MaleDDID[grepl("ASC",autism_counts$study)]),
                              sum(autism_counts$N_FemaleDDID[grepl("ASC",autism_counts$study)]),
                              sum(autism_counts$N_FemaleDDID[grepl("ASC",autism_counts$study)]),
                              sum(autism_counts$N_MaleNonDDID[grepl("ASC",autism_counts$study)]),
                              sum(autism_counts$N_MaleNonDDID[grepl("ASC",autism_counts$study)]),
                              sum(autism_counts$N_FemaleNonDDID[grepl("ASC",autism_counts$study)]),
                              sum(autism_counts$N_FemaleNonDDID[grepl("ASC",autism_counts$study)]),
                              sum(autism_counts$N_MaleDDID[grepl("SPARK",autism_counts$study)]),
                              sum(autism_counts$N_MaleDDID[grepl("SPARK",autism_counts$study)]),
                              sum(autism_counts$N_FemaleDDID[grepl("SPARK",autism_counts$study)]),
                              sum(autism_counts$N_FemaleDDID[grepl("SPARK",autism_counts$study)]),
                              sum(autism_counts$N_MaleNonDDID[grepl("SPARK",autism_counts$study)]),
                              sum(autism_counts$N_MaleNonDDID[grepl("SPARK",autism_counts$study)]),
                              sum(autism_counts$N_FemaleNonDDID[grepl("SPARK",autism_counts$study)]),
                              sum(autism_counts$N_FemaleNonDDID[grepl("SPARK",autism_counts$study)]),
                              sum(autism_counts$N_MaleDDID[grepl("GeneDx",autism_counts$study)]),
                              sum(autism_counts$N_MaleDDID[grepl("GeneDx",autism_counts$study)]),
                              sum(autism_counts$N_FemaleDDID[grepl("GeneDx",autism_counts$study)]),
                              sum(autism_counts$N_FemaleDDID[grepl("GeneDx",autism_counts$study)]),
                              sum(autism_counts$N_MaleNonDDID[grepl("GeneDx",autism_counts$study)]),
                              sum(autism_counts$N_MaleNonDDID[grepl("GeneDx",autism_counts$study)]),
                              sum(autism_counts$N_FemaleNonDDID[grepl("GeneDx",autism_counts$study)]),
                              sum(autism_counts$N_FemaleNonDDID[grepl("GeneDx",autism_counts$study)]),
                              sum(autism_counts$N_MaleDDID),
                              sum(autism_counts$N_MaleNonDDID),
                              sum(autism_counts$N_FemaleDDID),
                              sum(autism_counts$N_FemaleNonDDID),
                              sum(autism_counts$N_MaleDDID),
                              sum(autism_counts$N_MaleNonDDID),
                              sum(autism_counts$N_FemaleDDID),
                              sum(autism_counts$N_FemaleNonDDID))



autism_loop_vars$subsets = list(c("PTV_Proband_Male_DDID1","PTV_Proband_Female_DDID1","PTV_Proband_Male_DDID0","PTV_Proband_Female_DDID0"),
                                c("PTV_Sibling_Male_DDID1","PTV_Sibling_Female_DDID1","PTV_Sibling_Male_DDID0","PTV_Sibling_Female_DDID0"),
                                c("PTV_Proband_Male_DDID1","PTV_Proband_Female_DDID1","PTV_Proband_Male_DDID0","PTV_Proband_Female_DDID0"),
                                c("PTV_Sibling_Male_DDID1","PTV_Sibling_Female_DDID1","PTV_Sibling_Male_DDID0","PTV_Sibling_Female_DDID0"),
                                c("PTV_Proband_Male_DDID1","PTV_Proband_Female_DDID1","PTV_Proband_Male_DDID0","PTV_Proband_Female_DDID0"),
                                c("PTV_Sibling_Male_DDID1","PTV_Sibling_Female_DDID1","PTV_Sibling_Male_DDID0","PTV_Sibling_Female_DDID0"),
                                c("Mis2_Proband_Male_DDID1","Mis2_Proband_Female_DDID1","Mis2_Proband_Male_DDID0","Mis2_Proband_Female_DDID0"),
                                c("Mis2_Sibling_Male_DDID1","Mis2_Sibling_Female_DDID1","Mis2_Sibling_Male_DDID0","Mis2_Sibling_Female_DDID0"),
                                c("Mis2_Proband_Male_DDID1","Mis2_Proband_Female_DDID1","Mis2_Proband_Male_DDID0","Mis2_Proband_Female_DDID0"),
                                c("Mis2_Sibling_Male_DDID1","Mis2_Sibling_Female_DDID1","Mis2_Sibling_Male_DDID0","Mis2_Sibling_Female_DDID0"),
                                c("Mis2_Proband_Male_DDID1","Mis2_Proband_Female_DDID1","Mis2_Proband_Male_DDID0","Mis2_Proband_Female_DDID0"),
                                c("Mis2_Sibling_Male_DDID1","Mis2_Sibling_Female_DDID1","Mis2_Sibling_Male_DDID0","Mis2_Sibling_Female_DDID0"),
                                c("Mis1_Proband_Male_DDID1","Mis1_Proband_Female_DDID1","Mis1_Proband_Male_DDID0","Mis1_Proband_Female_DDID0"),
                                c("Mis1_Sibling_Male_DDID1","Mis1_Sibling_Female_DDID1","Mis1_Sibling_Male_DDID0","Mis1_Sibling_Female_DDID0"),
                                c("Mis1_Proband_Male_DDID1","Mis1_Proband_Female_DDID1","Mis1_Proband_Male_DDID0","Mis1_Proband_Female_DDID0"),
                                c("Mis1_Sibling_Male_DDID1","Mis1_Sibling_Female_DDID1","Mis1_Sibling_Male_DDID0","Mis1_Sibling_Female_DDID0"),
                                c("Mis1_Proband_Male_DDID1","Mis1_Proband_Female_DDID1","Mis1_Proband_Male_DDID0","Mis1_Proband_Female_DDID0"),
                                c("Mis1_Sibling_Male_DDID1","Mis1_Sibling_Female_DDID1","Mis1_Sibling_Male_DDID0","Mis1_Sibling_Female_DDID0"),
                                c("Mis0_Proband_Male_DDID1","Mis0_Proband_Female_DDID1","Mis0_Proband_Male_DDID0","Mis0_Proband_Female_DDID0"),
                                c("Mis0_Sibling_Male_DDID1","Mis0_Sibling_Female_DDID1","Mis0_Sibling_Male_DDID0","Mis0_Sibling_Female_DDID0"),
                                c("Mis0_Proband_Male_DDID1","Mis0_Proband_Female_DDID1","Mis0_Proband_Male_DDID0","Mis0_Proband_Female_DDID0"),
                                c("Mis0_Sibling_Male_DDID1","Mis0_Sibling_Female_DDID1","Mis0_Sibling_Male_DDID0","Mis0_Sibling_Female_DDID0"),
                                c("Mis0_Proband_Male_DDID1","Mis0_Proband_Female_DDID1","Mis0_Proband_Male_DDID0","Mis0_Proband_Female_DDID0"),
                                c("Mis0_Sibling_Male_DDID1","Mis0_Sibling_Female_DDID1","Mis0_Sibling_Male_DDID0","Mis0_Sibling_Female_DDID0"),
                                c("Syn_Proband_Male_DDID1","Syn_Proband_Female_DDID1","Syn_Proband_Male_DDID0","Syn_Proband_Female_DDID0"),
                                c("Syn_Sibling_Male_DDID1","Syn_Sibling_Female_DDID1","Syn_Sibling_Male_DDID0","Syn_Sibling_Female_DDID0"),
                                c("Syn_Proband_Male_DDID1","Syn_Proband_Female_DDID1","Syn_Proband_Male_DDID0","Syn_Proband_Female_DDID0"),
                                c("Syn_Sibling_Male_DDID1","Syn_Sibling_Female_DDID1","Syn_Sibling_Male_DDID0","Syn_Sibling_Female_DDID0"),
                                c("Syn_Proband_Male_DDID1","Syn_Proband_Female_DDID1","Syn_Proband_Male_DDID0","Syn_Proband_Female_DDID0"),
                                c("Syn_Sibling_Male_DDID1","Syn_Sibling_Female_DDID1","Syn_Sibling_Male_DDID0","Syn_Sibling_Female_DDID0"),
                                c("PTV_Proband_Male_DDID1","PTV_Proband_Female_DDID1","PTV_Proband_Male_DDID0","PTV_Proband_Female_DDID0"),
                                c("PTV_Proband_Male_DDID1","PTV_Proband_Male_DDID0"),
                                c("PTV_Proband_Female_DDID1","PTV_Proband_Female_DDID0"),
                                c("PTV_Proband_Male_DDID1","PTV_Proband_Female_DDID1","PTV_Proband_Male_DDID0","PTV_Proband_Female_DDID0"),
                                c("PTV_Proband_Male_DDID1","PTV_Proband_Male_DDID0"),
                                c("PTV_Proband_Female_DDID1","PTV_Proband_Female_DDID0"),
                                c("PTV_Proband_Male_DDID1","PTV_Proband_Female_DDID1","PTV_Proband_Male_DDID0","PTV_Proband_Female_DDID0"),
                                c("PTV_Proband_Male_DDID1","PTV_Proband_Male_DDID0"),
                                c("PTV_Proband_Female_DDID1","PTV_Proband_Female_DDID0"),
                                c("PTV_Proband_Male_DDID1","PTV_Proband_Male_DDID0"),
                                c("PTV_Proband_Female_DDID1","PTV_Proband_Female_DDID0"),
                                c("Mis2_Proband_Male_DDID1","Mis2_Proband_Male_DDID0"),
                                c("Mis2_Proband_Female_DDID1","Mis2_Proband_Female_DDID0"),
                                c("PTV_Proband_Male_DDID1","PTV_Proband_Female_DDID1"),
                                c("PTV_Proband_Male_DDID0","PTV_Proband_Female_DDID0"),
                                c("Mis2_Proband_Male_DDID1","Mis2_Proband_Female_DDID1"),
                                c("Mis2_Proband_Male_DDID0","Mis2_Proband_Female_DDID0"),
                                c("Mis1_Proband_Male_DDID1","Mis1_Proband_Female_DDID1"),
                                c("Mis1_Proband_Male_DDID0","Mis1_Proband_Female_DDID0"),
                                c("Mis0_Proband_Male_DDID1","Mis0_Proband_Female_DDID1"),
                                c("Mis0_Proband_Male_DDID0","Mis0_Proband_Female_DDID0"),
                                c("Syn_Proband_Male_DDID1","Syn_Proband_Female_DDID1"),
                                c("Syn_Proband_Male_DDID0","Syn_Proband_Female_DDID0"),
                                c("PTV_Proband_Male_DDID1","PTV_Proband_Female_DDID1","PTV_Proband_Male_DDID0","PTV_Proband_Female_DDID0"),
                                c("Mis2_Proband_Male_DDID1","Mis2_Proband_Female_DDID1","Mis2_Proband_Male_DDID0","Mis2_Proband_Female_DDID0"),
                                c("Mis1_Proband_Male_DDID1","Mis1_Proband_Female_DDID1","Mis1_Proband_Male_DDID0","Mis1_Proband_Female_DDID0"),
                                c("Mis0_Proband_Male_DDID1","Mis0_Proband_Female_DDID1","Mis0_Proband_Male_DDID0","Mis0_Proband_Female_DDID0"),
                                c("Syn_Proband_Male_DDID1","Syn_Proband_Female_DDID1","Syn_Proband_Male_DDID0","Syn_Proband_Female_DDID0"),
                                c("PTV_Proband_Male_DDID1"),
                                c("Mis2_Proband_Male_DDID1"),
                                c("PTV_Proband_Female_DDID1"),
                                c("Mis2_Proband_Female_DDID1"),
                                c("PTV_Proband_Male_DDID0"),
                                c("Mis2_Proband_Male_DDID0"),
                                c("PTV_Proband_Female_DDID0"),
                                c("Mis2_Proband_Female_DDID0"),
                                c("PTV_Proband_Male_DDID1"),
                                c("Mis2_Proband_Male_DDID1"),
                                c("PTV_Proband_Female_DDID1"),
                                c("Mis2_Proband_Female_DDID1"),
                                c("PTV_Proband_Male_DDID0"),
                                c("Mis2_Proband_Male_DDID0"),
                                c("PTV_Proband_Female_DDID0"),
                                c("Mis2_Proband_Female_DDID0"),
                                c("PTV_Proband_Male_DDID1"),
                                c("Mis2_Proband_Male_DDID1"),
                                c("PTV_Proband_Female_DDID1"),
                                c("Mis2_Proband_Female_DDID1"),
                                c("PTV_Proband_Male_DDID0"),
                                c("Mis2_Proband_Male_DDID0"),
                                c("PTV_Proband_Female_DDID0"),
                                c("Mis2_Proband_Female_DDID0"),
                                c("PTV_Proband_Male_DDID1"),
                                c("PTV_Proband_Male_DDID0"),
                                c("PTV_Proband_Female_DDID1"),
                                c("PTV_Proband_Female_DDID0"),
                                c("Mis2_Proband_Male_DDID1"),
                                c("Mis2_Proband_Male_DDID0"),
                                c("Mis2_Proband_Female_DDID1"),
                                c("Mis2_Proband_Female_DDID0"))



autism_loop_vars$mutation_rate = c("mu_snp_PTV",
                                   "mu_snp_PTV",
                                   "mu_snp_PTV",
                                   "mu_snp_PTV",
                                   "mu_snp_PTV",
                                   "mu_snp_PTV",
                                   "mu_snp_Mis2",
                                   "mu_snp_Mis2",
                                   "mu_snp_Mis2",
                                   "mu_snp_Mis2",
                                   "mu_snp_Mis2",
                                   "mu_snp_Mis2",
                                   "mu_snp_Mis1",
                                   "mu_snp_Mis1",
                                   "mu_snp_Mis1",
                                   "mu_snp_Mis1",
                                   "mu_snp_Mis1",
                                   "mu_snp_Mis1",
                                   "mu_snp_Mis0",
                                   "mu_snp_Mis0",
                                   "mu_snp_Mis0",
                                   "mu_snp_Mis0",
                                   "mu_snp_Mis0",
                                   "mu_snp_Mis0",
                                   "mu_snp_Syn",
                                   "mu_snp_Syn",
                                   "mu_snp_Syn",
                                   "mu_snp_Syn",
                                   "mu_snp_Syn",
                                   "mu_snp_Syn",
                                   "mu_snp_PTV",
                                   "mu_snp_PTV",
                                   "mu_snp_PTV",
                                   "mu_snp_PTV",
                                   "mu_snp_PTV",
                                   "mu_snp_PTV",
                                   "mu_snp_PTV",
                                   "mu_snp_PTV",
                                   "mu_snp_PTV",
                                   "mu_snp_PTV",
                                   "mu_snp_PTV",
                                   "mu_snp_Mis2",
                                   "mu_snp_Mis2",
                                   "mu_snp_PTV",
                                   "mu_snp_PTV",
                                   "mu_snp_Mis2",
                                   "mu_snp_Mis2",
                                   "mu_snp_Mis1",
                                   "mu_snp_Mis1",
                                   "mu_snp_Mis0",
                                   "mu_snp_Mis0",
                                   "mu_snp_Syn",
                                   "mu_snp_Syn",
                                   "mu_snp_PTV",
                                   "mu_snp_Mis2",
                                   "mu_snp_Mis1",
                                   "mu_snp_Mis0",
                                   "mu_snp_Syn",
                                   "mu_snp_PTV",
                                   "mu_snp_Mis2",
                                   "mu_snp_PTV",
                                   "mu_snp_Mis2",
                                   "mu_snp_PTV",
                                   "mu_snp_Mis2",
                                   "mu_snp_PTV",
                                   "mu_snp_Mis2",
                                   "mu_snp_PTV",
                                   "mu_snp_Mis2",
                                   "mu_snp_PTV",
                                   "mu_snp_Mis2",
                                   "mu_snp_PTV",
                                   "mu_snp_Mis2",
                                   "mu_snp_PTV",
                                   "mu_snp_Mis2",
                                   "mu_snp_PTV",
                                   "mu_snp_Mis2",
                                   "mu_snp_PTV",
                                   "mu_snp_Mis2",
                                   "mu_snp_PTV",
                                   "mu_snp_Mis2",
                                   "mu_snp_PTV",
                                   "mu_snp_Mis2",
                                   "mu_snp_PTV",
                                   "mu_snp_PTV",
                                   "mu_snp_PTV",
                                   "mu_snp_PTV",
                                   "mu_snp_Mis2",
                                   "mu_snp_Mis2",
                                   "mu_snp_Mis2",
                                   "mu_snp_Mis2"
)




baseprev = 0.0276

autism_loop_vars$prevalences = c(rep(baseprev, 30),
                                 baseprev,
                                 baseprev * (0.0430/0.0276),
                                 baseprev * (0.0114/0.0276),
                                 baseprev,
                                 baseprev * (0.0430/0.0276),
                                 baseprev * (0.0114/0.0276),
                                 baseprev,
                                 baseprev * (0.0430/0.0276),
                                 baseprev * (0.0114/0.0276),
                                 baseprev * (0.0430/0.0276),
                                 baseprev * (0.0114/0.0276),
                                 baseprev * (0.0430/0.0276),
                                 baseprev * (0.0114/0.0276),
                                 baseprev * 0.379,
                                 baseprev * (1 - 0.379),
                                 baseprev * 0.379,
                                 baseprev * (1 - 0.379),
                                 baseprev * 0.379,
                                 baseprev * (1 - 0.379),
                                 baseprev * 0.379,
                                 baseprev * (1 - 0.379),
                                 baseprev * 0.379,
                                 baseprev * (1 - 0.379),
                                 baseprev,
                                 baseprev,
                                 baseprev,
                                 baseprev,
                                 baseprev,
                                 baseprev * (0.0430/0.0276)* 0.369,
                                 baseprev * (0.0430/0.0276)* 0.369,
                                 baseprev * (0.0114/0.0276)* 0.421,
                                 baseprev * (0.0114/0.0276)* 0.421,
                                 baseprev * (0.0430/0.0276)* (1 - 0.369),
                                 baseprev * (0.0430/0.0276)* (1 - 0.369),
                                 baseprev * (0.0114/0.0276)* (1 - 0.421),
                                 baseprev * (0.0114/0.0276)* (1 - 0.421),
                                 baseprev * (0.0430/0.0276)* 0.369,
                                 baseprev * (0.0430/0.0276)* 0.369,
                                 baseprev * (0.0114/0.0276)* 0.421,
                                 baseprev * (0.0114/0.0276)* 0.421,
                                 baseprev * (0.0430/0.0276)* (1 - 0.369),
                                 baseprev * (0.0430/0.0276)* (1 - 0.369),
                                 baseprev * (0.0114/0.0276)* (1 - 0.421),
                                 baseprev * (0.0114/0.0276)* (1 - 0.421),
                                 baseprev * (0.0430/0.0276)* 0.369,
                                 baseprev * (0.0430/0.0276)* 0.369,
                                 baseprev * (0.0114/0.0276)* 0.421,
                                 baseprev * (0.0114/0.0276)* 0.421,
                                 baseprev * (0.0430/0.0276)* (1 - 0.369),
                                 baseprev * (0.0430/0.0276)* (1 - 0.369),
                                 baseprev * (0.0114/0.0276)* (1 - 0.421),
                                 baseprev * (0.0114/0.0276)* (1 - 0.421),
                                 baseprev * (0.0430/0.0276)* 0.369,
                                 baseprev * (0.0430/0.0276)* (1 - 0.369),
                                 baseprev * (0.0114/0.0276)* 0.421,
                                 baseprev * (0.0114/0.0276)* (1 - 0.421),
                                 baseprev * (0.0430/0.0276)* 0.369,
                                 baseprev * (0.0430/0.0276)* (1 - 0.369),
                                 baseprev * (0.0114/0.0276)* 0.421,
                                 baseprev * (0.0114/0.0276)* (1 - 0.421)


)

# Remove analyses that are not used by the final visualization workflow.
unused_model_names <- c(
  grep("^(New|Fu2022)", autism_loop_vars$names, value = TRUE),
  "Combined DDID Mis1", "Combined Non-DDID Mis1",
  "Combined DDID Mis0", "Combined Non-DDID Mis0",
  "Combined DDID Syn", "Combined Non-DDID Syn",
  "SPARK + ASC PTV", "SPARK + ASC Mis2", "SPARK + ASC Mis1",
  "SPARK + ASC Mis0", "SPARK + ASC Syn"
)

loop_var_lengths <- vapply(autism_loop_vars, length, integer(1))
if (length(unique(loop_var_lengths)) != 1) {
  stop("Autism loop-variable fields have different lengths: ",
       paste(names(loop_var_lengths), loop_var_lengths, sep = "=", collapse = ", "))
}

keep_model <- !autism_loop_vars$names %in% unused_model_names
autism_loop_vars <- lapply(autism_loop_vars, function(x) x[keep_model])

# These entries define useful count contrasts for get_genetic_data(), but the
# visualization rescales the combined models rather than using separate fits.
count_contrast_names <- c(
  "Combined ASC PTV", "Male ASC PTV", "Female ASC PTV",
  "Combined SPARK PTV", "Male SPARK PTV", "Female SPARK PTV",
  "Combined GeneDx PTV", "Male GeneDx PTV", "Female GeneDx PTV"
)
autism_loop_vars$fit_model <- !autism_loop_vars$names %in% count_contrast_names


# Derive prevalence weights from the current pedigree data rather than stale
# sample-size constants.
N_SPARK = sum(autism_counts$N_Proband[grepl("SPARK", autism_counts$study)])
N_ASC = sum(autism_counts$N_Proband[grepl("ASC", autism_counts$study)])
N_GeneDX = sum(autism_counts$N_Proband[grepl("GeneDx", autism_counts$study)])
N_DDD = sum(autism_counts$N_Proband[grepl("DDD", autism_counts$study)])
N_total = N_SPARK + N_ASC + N_GeneDX + N_DDD

weight_SPARK = N_SPARK/N_total
weight_ASC = N_ASC/N_total
weight_GeneDX = N_GeneDX/N_total
weight_DDD = N_DDD/N_total

prev_weighted = (weight_SPARK * 0.0276) + (weight_ASC * 0.0276) + (weight_GeneDX * 0.01) + (weight_DDD * 0.01)

prev_factors = c(prev_weighted/0.0276,0.01/0.0276,1)


autism_data <- list(counts = autism_counts,
                    features = autism_features,
                    loop_vars = autism_loop_vars,
                    ptv_scale_factor = scale_ptv_factor,
                    baseprev = baseprev,
                    prev_factors = prev_factors)
