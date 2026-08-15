mutvar_enrichment_table <- function(data,output) {
  mutvar_est_df <- data.frame()
  enrich_est_df <- data.frame()

  for (prevfactornumber in 1:length(data$prev_factors)) {

    prev_factor = data$prev_factors[prevfactornumber]
    print(prev_factor)
    for (outputnum in 1:length(data$loop_vars$names)) {

      prev_nonmod = data$loop_vars$prevalences[outputnum]
      prev_mod = data$loop_vars$prevalences[outputnum] * prev_factor

      variant_class = if (grepl("PTV",data$loop_vars$subsets[outputnum])) {
        "PTV"
      } else if (grepl("Mis2",data$loop_vars$subsets[outputnum])) {
        "Mis2"
      }else if (grepl("Mis1",data$loop_vars$subsets[outputnum])) {
        "Mis1"
      }else if (grepl("Mis0",data$loop_vars$subsets[outputnum])) {
        "Mis0"
      }else if (grepl("Syn",data$loop_vars$subsets[outputnum])) {
        "Syn"
      }
      role = if (grepl("Proband",data$loop_vars$subsets[outputnum])) {
        "Proband"
      } else if (grepl("Sibling",data$loop_vars$subsets[outputnum])) {
        "Sibling"
      }

      if (all(is.na(output[[prevfactornumber]][[outputnum]]))) {
        next
      }

      model = output[[prevfactornumber]][[outputnum]]
      #model = BurdenMLE_DN_models[[outputnum]]

      if (variant_class == "PTV") {
        mutvar = model$mutvar_output$total_mutvar * data$ptv_scale_factor
        print(mutvar)
        mutvar_lower = model$mutvar_output$mutvar_CI[1] * data$ptv_scale_factor
        mutvar_upper = model$mutvar_output$mutvar_CI[2] * data$ptv_scale_factor

      } else {
        mutvar = model$mutvar_output$total_mutvar
        mutvar_lower = model$mutvar_output$mutvar_CI[1]
        mutvar_upper = model$mutvar_output$mutvar_CI[2]
      }
      mutvar_p = model$mutvar_output$mutvar_p
      if (is.null(mutvar_p)) {mutvar_p <- NA}

      peneff = model$penetrance$effective_penetrance
      peneff_lower = model$penetrance$effective_penetrance_CI[1]
      peneff_upper = model$penetrance$effective_penetrance_CI[2]


      output_mutvar_df <- data.frame(name = data$loop_vars$names[outputnum],
                                 variant_class = variant_class,
                                 prev_nonmod = prev_nonmod,
                                 prev_mod = prev_mod,
                                 prev_factor = prev_factor,
                                 role = role,
                                 mutvar = mutvar,
                                 mutvar_lower = mutvar_lower,
                                 mutvar_upper = mutvar_upper,
                                 mutvar_p = mutvar_p,
                                 peneff = peneff,
                                 peneff_lower = peneff_lower,
                                 peneff_upper = peneff_upper)
      mutvar_est_df = rbind(mutvar_est_df,output_mutvar_df)

      mutvar_enrich = model$mutvar_output$enrichment
      mutvar_lower = model$mutvar_output$enrich_CI[1,]
      mutvar_upper = model$mutvar_output$enrich_CI[2,]

      frac_mutvar = model$mutvar_output$frac_mutvar
      frac_expected = model$mutvar_output$frac_expected

      output_enrich_df <- data.frame(name = data$loop_vars$names[outputnum],
                                     variant_class = variant_class,
                                     prev_nonmod = prev_nonmod,
                                     prev_mod = prev_mod,
                                     prev_factor = prev_factor,
                                     role = role,
                                     annot = c("LOEUF1_mu1","LOEUF1_mu2",paste0("LOEUF",2:5)),
                                     mutvar_enrich = mutvar_enrich,
                                     mutvar_enrich_lower = mutvar_lower,
                                     mutvar_enrich_upper = mutvar_upper,
                                     frac_mutvar = frac_mutvar,
                                     frac_expected = frac_expected)

      enrich_est_df = rbind(enrich_est_df, output_enrich_df)
    }

  }
  return(list(mutvar = mutvar_est_df,
              enrichment = enrich_est_df))
}

bootstrap_function <- function(model, function_to_bootstrap, ...) {
  args <- list(...)  # Capture additional arguments as a list

  lapply(seq_along(model$bootstrap_output$bootstrap_delta), function(iter) {
    # Start every replicate from the original arguments. Otherwise bootstrap
    # indices after the first replicate are applied to already-resampled data.
    iter_args <- args
    model_boot <- model
    model_boot$conditional_likelihood <- model_boot$conditional_likelihood[
      model$bootstrap_output$bootstrap_indices[, iter], , drop = FALSE
    ]
    model_boot$features <- model_boot$features[
      model$bootstrap_output$bootstrap_indices[, iter], , drop = FALSE
    ]
    model_boot$delta <- model$bootstrap_output$bootstrap_delta[[iter]]

    # If genetic_data is provided in ..., subset it
    if ("genetic_data" %in% names(iter_args)) {
      genetic_data_boot <- iter_args$genetic_data[
        model$bootstrap_output$bootstrap_indices[, iter], , drop = FALSE
      ]
      iter_args$genetic_data <- genetic_data_boot
    }

    # Call the function with the modified arguments
    do.call(function_to_bootstrap, c(list(model = model_boot), iter_args))
  })
}

get_fraccase <- function(model,
                         genetic_data,
                         RR_thresh,
                         gamma_scaling_factor = 1) {
  sum(2*genetic_data$case_rate*posterior_expectation(model,
                                                     genetic_data,
                                                     function(x) {
                                                       exp(gamma_scaling_factor*x) * (exp(x) > RR_thresh)
                                                     },
                                                     grid_size = 10))
}

get_fraccase_df <- function(data,
                            modelPTV,genetic_dataPTV,
                            modelMis2, genetic_dataMis2,
                            gamma_scaling_factor = 1,
                            boot = TRUE,
                            RR_range = c(2:20)) {
  fraccases_effsizethresh_df = data.frame()

  for (RR_thresh in RR_range) {
    print(RR_thresh)
    frac_casesgreater_pergene_PTV =  get_fraccase(model = modelPTV,
                                                  genetic_data = genetic_dataPTV,
                                                  gamma_scaling_factor = gamma_scaling_factor,
                                                  RR_thresh = RR_thresh)

    frac_casesgreater_pergene_Mis2 =  get_fraccase(model = modelMis2,
                                                   genetic_data = genetic_dataMis2,
                                                   gamma_scaling_factor = gamma_scaling_factor,
                                                   RR_thresh = RR_thresh)
    frac_casesgreater_combined = frac_casesgreater_pergene_PTV*data$ptv_scale_factor + frac_casesgreater_pergene_Mis2


    #bootstrap
    if (boot) {
      PTV_boot <- unlist(bootstrap_function(modelPTV,
                                            get_fraccase,
                                            genetic_data = genetic_dataPTV,
                                            RR_thresh = RR_thresh,
                                            gamma_scaling_factor = gamma_scaling_factor))
      Mis2_boot <- unlist(bootstrap_function(modelMis2,
                                             get_fraccase,
                                             genetic_data = genetic_dataMis2,
                                             RR_thresh = RR_thresh,
                                             gamma_scaling_factor = gamma_scaling_factor))
      frac_casesgreater_combined_boot = PTV_boot*data$ptv_scale_factor + Mis2_boot
    } else {
      PTV_boot <- rep(NA, length(modelPTV$bootstrap_output$bootstrap_delta))
      Mis2_boot <- rep(NA, length(modelPTV$bootstrap_output$bootstrap_delta))
      frac_casesgreater_combined_boot = PTV_boot*data$ptv_scale_factor + Mis2_boot
    }

    iter_df <- data.frame(RR_thresh = RR_thresh,
                          fraccases_greater_PTV = frac_casesgreater_pergene_PTV,
                          fraccases_greater_PTV_lowerCI = quantile(PTV_boot,0.025,na.rm = TRUE),
                          fraccases_greater_PTV_upperCI = quantile(PTV_boot,0.975,na.rm = TRUE),
                          fraccases_greater_Mis2 = frac_casesgreater_pergene_Mis2,
                          fraccases_greater_Mis2_lowerCI = quantile(Mis2_boot,0.025,na.rm = TRUE),
                          fraccases_greater_Mis2_upperCI = quantile(Mis2_boot,0.975,na.rm = TRUE),
                          fraccases_greater_PTVplusMis2 = frac_casesgreater_combined,
                          fraccases_greater_PTVplusMis2_lowerCI = quantile(frac_casesgreater_combined_boot, 0.025,na.rm = TRUE),
                          fraccases_greater_PTVplusMis2_upperCI = quantile(frac_casesgreater_combined_boot, 0.975,na.rm = TRUE),
                          gamma_scaling_factor = gamma_scaling_factor)

    fraccases_effsizethresh_df = rbind(fraccases_effsizethresh_df,iter_df)

  }
  return(fraccases_effsizethresh_df)
}

get_fraccase_df_general <- function(variant_classes,
                                    gamma_scaling_factor = 1,
                                    boot = TRUE,
                                    RR_range = 2:20) {
  fraccases_effsizethresh_df = data.frame()

  for (RR_thresh in RR_range) {
    print(RR_thresh)

    results <- lapply(names(variant_classes), function(name) {
      vc <- variant_classes[[name]]

      fraccase <- get_fraccase(model = vc$model,
                               genetic_data = vc$genetic_data,
                               gamma_scaling_factor = gamma_scaling_factor,
                               RR_thresh = RR_thresh)

      if (boot) {
        boot_vals <- unlist(bootstrap_function(vc$model,
                                               get_fraccase,
                                               genetic_data = vc$genetic_data,
                                               RR_thresh = RR_thresh,
                                               gamma_scaling_factor = gamma_scaling_factor))
      } else {
        boot_vals <- rep(NA, length(vc$model$bootstrap_output$bootstrap_delta))
      }

      list(name = name,
           fraccase = fraccase,
           boot = boot_vals,
           scale = ifelse(is.null(vc$scale), 1, vc$scale))
    })

    # Combine results across classes
    combined_fraccase <- Reduce(`+`, lapply(results, function(r) r$fraccase * r$scale))
    combined_boot <- Reduce(`+`, lapply(results, function(r) r$boot * r$scale))

    # Create data.frame row
    iter_row <- list(RR_thresh = RR_thresh, gamma_scaling_factor = gamma_scaling_factor)

    for (res in results) {
      iter_row[[paste0("fraccases_greater_", res$name)]] <- res$fraccase
      iter_row[[paste0("fraccases_greater_", res$name, "_lowerCI")]] <- quantile(res$boot, 0.025, na.rm = TRUE)
      iter_row[[paste0("fraccases_greater_", res$name, "_upperCI")]] <- quantile(res$boot, 0.975, na.rm = TRUE)
    }

    iter_row[["fraccases_greater_combined"]] <- combined_fraccase
    iter_row[["fraccases_greater_combined_lowerCI"]] <- quantile(combined_boot, 0.025, na.rm = TRUE)
    iter_row[["fraccases_greater_combined_upperCI"]] <- quantile(combined_boot, 0.975, na.rm = TRUE)

    fraccases_effsizethresh_df <- rbind(fraccases_effsizethresh_df, as.data.frame(iter_row))
  }

  return(fraccases_effsizethresh_df)
}



ngenes_for_fraccases <- function(data,
                                 modelPTV,genetic_dataPTV,
                                 modelMis2, genetic_dataMis2,
                                 gamma_scaling_factor = 1) {

  ngenes_for_fraccases <- sapply(2:20,
                                 function(RR_thresh) {
                                   frac_casesgreater_pergene_PTV =  (2*genetic_dataPTV$case_rate*posterior_expectation(modelPTV,
                                                                                                                       genetic_dataPTV,
                                                                                                                       function(x) {
                                                                                                                         exp(gamma_scaling_factor*x) * (exp(x) > RR_thresh)
                                                                                                                       },
                                                                                                                       grid_size = 10))

                                   frac_casesgreater_pergene_Mis2 =  (2*genetic_dataMis2$case_rate*posterior_expectation(modelMis2,
                                                                                                                         genetic_dataMis2,
                                                                                                                         function(x) {
                                                                                                                           exp(gamma_scaling_factor*x) * (exp(x) > RR_thresh)
                                                                                                                         },
                                                                                                                         grid_size = 10))

                                   fraccase_greater_combined = (frac_casesgreater_pergene_PTV*scale_ptv_factor) + frac_casesgreater_pergene_Mis2

                                   fraccase_greater_combined_ordered = fraccase_greater_combined[order(fraccase_greater_combined, decreasing = TRUE)]
                                   return(cumsum(fraccase_greater_combined_ordered))
                                 })

  colnames(ngenes_for_fraccases) = paste0("RR>",2:20)
  rownames(ngenes_for_fraccases) <- 1:nrow(ngenes_for_fraccases)
  return(ngenes_for_fraccases)
}

countgene_func <- function(model,
                           genetic_data,
                           RR_thresh,
                           filter) {
  sum(posterior_expectation(model,
                            genetic_data,
                            function(x) {exp(x) > RR_thresh},
                            grid_size = 10)[filter])
}

count_genes_by_effectsize <- function(model, genetic_data, genes_to_analyze = NULL, effsize_range =  c(2,5,10,20), boot = FALSE) {

  if (!is.null(genes_to_analyze) & boot) {
    stop("Cannot Subset and Bootstrap at the Same Time")
  }
  if (is.null(genes_to_analyze)) {
    filter = rep(TRUE, nrow(genetic_data))
  } else {
    filter = rownames(genetic_data) %in% genes_to_analyze
  }
  effsizethresh_df = data.frame()

  for (RR_thresh in effsize_range) {
    print(RR_thresh)
    num_genes =  countgene_func(model,
                                genetic_data,
                                RR_thresh,
                                filter)


    if (boot) {
      num_genes_boot <- unlist(bootstrap_function(model,
                                            countgene_func,
                                            genetic_data = genetic_data,
                                            RR_thresh = RR_thresh))
      num_genes_CI = quantile(num_genes_boot, c(0.025, 0.975))
      iter_df <- data.frame(RR_thresh = RR_thresh,
                            num_genes = num_genes,
                            num_genes_lower = num_genes_CI[1],
                            num_genes_upper = num_genes_CI[2])
    } else {
      iter_df <- data.frame(RR_thresh = RR_thresh,
                            num_genes = num_genes)
    }


    effsizethresh_df = rbind(effsizethresh_df, iter_df)

  }
  return(effsizethresh_df)
}


mutvar_per_gene <- function(model,genetic_data,prevalence,gamma_scaling_factor) {
  ((prevalence * 2 * genetic_data$case_rate)/(1 - prevalence)) * posterior_expectation(model,
                                                                                       genetic_data,
                                                                                       function(x) {
                                                                                         (exp(gamma_scaling_factor*x) -1)^2
                                                                                       },
                                                                                       grid_size = 10)
}

aggregate_mutvar_per_gene <- function(data,
                                      modelPTV,genetic_dataPTV,
                                      modelMis2, genetic_dataMis2,
                                      prevalence,
                                      gamma_scaling_factor = 1) {


  mutvar_pergene_PTV <- mutvar_per_gene(modelPTV,genetic_dataPTV,prevalence,gamma_scaling_factor)

  mutvar_pergene_Mis2 <- mutvar_per_gene(modelMis2,genetic_dataMis2,prevalence,gamma_scaling_factor)


  mutvar_pergene_combined = (mutvar_pergene_PTV*data$ptv_scale_factor) + mutvar_pergene_Mis2
  names(mutvar_pergene_combined) <- rownames(genetic_dataPTV)

  mutvar_pergene_ordered = mutvar_pergene_combined[order(mutvar_pergene_combined, decreasing = TRUE)]

  #bootstrap
  mutvar_pergene_PTV_boot <- bootstrap_function(modelPTV,
                                                mutvar_per_gene,
                                                genetic_data = genetic_dataPTV,
                                                prevalence = prevalence,
                                                gamma_scaling_factor = 1)

  mutvar_pergene_Mis2_boot <- bootstrap_function(modelMis2,
                                                 mutvar_per_gene,
                                                 genetic_data = genetic_dataMis2,
                                                 prevalence = prevalence,
                                                 gamma_scaling_factor = 1)

  mutvar_pergene_combined_boot <- sapply(1:length(mutvar_pergene_PTV_boot),
                                         function(i) {
                                           mutvar_pergene_combined_boot = (mutvar_pergene_PTV_boot[[i]]*data$ptv_scale_factor) + mutvar_pergene_Mis2_boot[[i]]
                                           names(mutvar_pergene_combined_boot) <- rownames(genetic_dataPTV)

                                           mutvar_pergene_ordered_boot = mutvar_pergene_combined_boot[order(mutvar_pergene_combined_boot, decreasing = TRUE)]
                                         })
  rownames(mutvar_pergene_combined_boot) <- 1:nrow(mutvar_pergene_combined_boot)
  return(list(full = mutvar_pergene_ordered,
              boot = mutvar_pergene_combined_boot))
}

BurdenMLE_DN_power_forecasting <- function(old_genetic_data,
                                       model,
                                       N_new_case,
                                       gamma_scaling_factor = 1,
                                       effect_size_thresholds = c(0,2,5,10,20)) {

  #sample effect sizes from posteriors
  log_gamma_samples = posterior_gene_samples(model) * gamma_scaling_factor

  #Sample some new counts
  new_case_counts <- rpois(nrow(old_genetic_data), 2 * N_new_case * exp(log_gamma_samples) * old_genetic_data$case_rate)

  #combine
  combined_counts = new_case_counts + old_genetic_data$case_count

  new_genetic_data = data.frame(case_count = combined_counts,
                                case_rate = old_genetic_data$case_rate,
                                N = old_genetic_data$N +N_new_case)
  rownames(new_genetic_data) = rownames(old_genetic_data)
  new_genetic_data$expected_count = 2 * new_genetic_data$N * new_genetic_data$case_rate

  #Estimate some poisson P values
  poisson_p <- ppois(combined_counts-0.001,new_genetic_data$expected_count, lower.tail = FALSE)

  #Sample some null counts
  null_counts = rpois(nrow(new_genetic_data), new_genetic_data$expected_count)

  #Get some null P values
  null_p <- ppois(null_counts-0.001,new_genetic_data$expected_count, lower.tail = FALSE)


  # QQ plot
  qq <- ggplot(mapping = aes(y = -log10(poisson_p[order(poisson_p,decreasing = TRUE)]),
                             x = -log10(null_p[order(null_p,decreasing = TRUE)])))+
    geom_point()+
    geom_abline()+
    geom_hline(yintercept = -log10( 0.05/length(poisson_p)), linetype = "dashed")+
    labs(x = "Simulated Null Quantiles", y = "Simulated Effect Quantiles")+
    theme_bw()

  #How many significant genes?
  count_table = data.frame()
  for (RR_thresh in effect_size_thresholds) {

    gene_greaterthan_thresh =posterior_expectation(model,
                                                   old_genetic_data,
                                                   function(x) {exp(x) > RR_thresh},
                                                   grid_size = 10)

    iter_df <- data.frame(threshold = RR_thresh,
                          count_fdr = sum((p.adjust(poisson_p, method = "fdr") < 0.05) * gene_greaterthan_thresh),
                          count_bonferroni = sum((poisson_p < 0.05/length(poisson_p)) * gene_greaterthan_thresh))
    count_table <- rbind(count_table,iter_df)
  }

  return(list(count_table = count_table,
              qq = qq))

}

BurdenMLE_DN_power_forecasting_mod <- function(old_genetic_data,
                                           model,
                                           N_new_case,
                                           gamma_scaling_factor = 1,
                                           effect_size_thresholds = c(0,2,5,10,20)) {

  #sample effect sizes from posteriors. Define these to be the "true effect sizes"
  log_gamma_samples = posterior_gene_samples(model)

  #Scale the effect sizes for the new counts
  log_gamma_samples_scale =log_gamma_samples * gamma_scaling_factor

  #Sample some new counts
  new_case_counts <- rpois(nrow(old_genetic_data), 2 * N_new_case * exp(log_gamma_samples_scale) * old_genetic_data$case_rate)

  #combine
  combined_counts = new_case_counts + old_genetic_data$case_count

  new_genetic_data = data.frame(case_count = combined_counts,
                                case_rate = old_genetic_data$case_rate,
                                N = old_genetic_data$N +N_new_case)
  rownames(new_genetic_data) = rownames(old_genetic_data)
  new_genetic_data$expected_count = 2 * new_genetic_data$N * new_genetic_data$case_rate

  #Estimate some poisson P values
  poisson_p <- ppois(combined_counts-0.001,new_genetic_data$expected_count, lower.tail = FALSE)

  #Sample some null counts
  null_counts = rpois(nrow(new_genetic_data), new_genetic_data$expected_count)

  #Get some null P values
  null_p <- ppois(null_counts-0.001,new_genetic_data$expected_count, lower.tail = FALSE)


  # QQ plot
  qq <- ggplot(mapping = aes(y = -log10(poisson_p[order(poisson_p,decreasing = TRUE)]),
                             x = -log10(null_p[order(null_p,decreasing = TRUE)])))+
    geom_point()+
    geom_abline()+
    geom_hline(yintercept = -log10( 0.05/length(poisson_p)), linetype = "dashed")+
    labs(x = "Simulated Null Quantiles", y = "Simulated Effect Quantiles")+
    theme_bw()

  #How many significant genes?
  count_table = data.frame()
  for (RR_thresh in effect_size_thresholds) {

    #Rather than using the posterior distribution, simply use the sampled effect size as the 'true' effect size

    gene_greaterthan_thresh = exp(log_gamma_samples) >= RR_thresh
    iter_df <- data.frame(threshold = RR_thresh,
                          count_fdr = sum((p.adjust(poisson_p, method = "fdr") < 0.05) * gene_greaterthan_thresh),
                          count_bonferroni = sum((poisson_p < 0.05/length(poisson_p)) * gene_greaterthan_thresh))
    count_table <- rbind(count_table,iter_df)
  }

  return(list(count_table = count_table,
              qq = qq))

}

make_supptable <- function(tables,names_to_keep) {
  mutvar_table = tables$mutvar
  enrich_table = tables$enrichment

  mutvar_table = mutvar_table[mutvar_table$name %in% names_to_keep,]
  enrich_table = enrich_table[enrich_table$name %in% names_to_keep,]

  # Reshape enrich_table to wide format
  enrich_wide <- reshape(enrich_table,
                         idvar = c("name", "variant_class", "prev_nonmod", "prev_mod", "prev_factor", "role"),
                         timevar = "annot",
                         direction = "wide")

  # Rename columns to follow the new naming convention
  colnames(enrich_wide) <- gsub("mutvar_enrich\\.", "Enrichment_", colnames(enrich_wide))
  colnames(enrich_wide) <- gsub("mutvar_enrich_lower\\.", "Enrichment_Lower95CI_", colnames(enrich_wide))
  colnames(enrich_wide) <- gsub("mutvar_enrich_upper\\.", "Enrichment_Upper95CI_", colnames(enrich_wide))
  colnames(enrich_wide) <- gsub("frac_mutvar\\.", "FractionMutVar_", colnames(enrich_wide))
  colnames(enrich_wide) <- gsub("frac_expected\\.", "FractionExpected_", colnames(enrich_wide))

  # Merge the reshaped enrich_wide with mutvar_table
  mutvar_table_merged <- merge(mutvar_table, enrich_wide,
                           by = c("name", "variant_class", "prev_nonmod", "prev_mod", "prev_factor", "role"),
                           all.x = TRUE)

  supptable = mutvar_table_merged[,c("name",
                                 "variant_class",
                                 "prev_mod",
                                 "role",
                                 "mutvar",
                                 "mutvar_lower",
                                 "mutvar_upper",
                                 "peneff",
                                 "peneff_lower",
                                 "peneff_upper",
                                 names(mutvar_table_merged)[16:43])]

  names(supptable) = c("Model Name",
                       "Variant Class",
                       "Prevalence",
                       "Role",
                       "MutVar",
                       "MutVar_lower95CI",
                       "MutVar_upper95CI",
                       "EffectivePenetrance",
                       "EffectivePenetrance_Lower95CI",
                       "EffectivePenetrance_Upper95CI",
                       names(mutvar_table_merged)[16:43])


  return(supptable)


}

#define function to get count tables
get_genetic_data <- function(i, data) {
  counts = data$counts
  subsets = data$loop_vars$subsets
  mutation_rate = data$loop_vars$mutation_rate
  N_subset = data$loop_vars$N_subset
  names = data$loop_vars$names
  print(names[i])

  case_count = rep(0, nrow(rowData(counts)))

  for (assay in subsets[[i]]) {
    if (length(data$loop_vars$datasets) > 0) {
      if (sum(data$loop_vars$datasets[[i]]) > 1) {
        assay_count = rowSums(assays(counts)[[assay]][,data$loop_vars$datasets[[i]]])
      } else {
        assay_count = assays(counts)[[assay]][,data$loop_vars$datasets[[i]]]
      }
    } else {
      assay_count = rowSums(assays(counts)[[assay]])

    }
    case_count = case_count + assay_count
  }

  mu = rowData(counts)[,mutation_rate[i]]
  #  mutation_rate = rowData(counts)$mut.ptv.gnomad.v2

  posterior_mu_factor = rowData(counts)$PosteriorMuCorrectionFactor

  N = N_subset[i]

  input_df <- data.frame(case_count = case_count,
                         case_rate = mu * posterior_mu_factor,
                         N = N)
  rownames(input_df) = rownames(rowData(counts))

  filter_impossiblecount = !(input_df$case_count > 0 & input_df$case_rate == 0)


  genetic_data = process_data_trio(input_df[filter_impossiblecount,],
                                   data$features[filter_impossiblecount,])

  return(list(genetic_data = genetic_data,
              features = data$features[filter_impossiblecount,]))

}

binomial_LL_withposterior <- function(c,data,totals,N_1,N_tot,model,genetic_data) {
  posterior_terms = posterior_expectation(model,genetic_data,function(x) {exp(x * (c-1))},10)
  probability_factors = (N_1/N_tot) * posterior_terms


  LLs = dbinom(data,totals,probability_factors,log = TRUE)
  return(sum(LLs))
}

binomial_analysis <- function(model,
                              genetic_data_total,
                              genetic_data_subsample,
                              param_grid = seq(0.6,1.35,length.out  = 30)) {

  counts_total = genetic_data_total$case_count
  counts_subsample = genetic_data_subsample$case_count
  N_total = genetic_data_total$N[1]
  N_subsample = genetic_data_subsample$N[1]
  print(N_subsample)

  LL_grid_posterior = sapply(param_grid,
                             binomial_LL_withposterior,
                             counts_subsample,
                             counts_total,
                             N_subsample,
                             N_total,
                             model,
                             genetic_data_total)

  return(list(LL_output = data.frame(val = param_grid,
                                     LL = LL_grid_posterior),
              MLE = param_grid[which.max(LL_grid_posterior)]))


}

get_scaled_mutvar <- function(model,
                                    genetic_data,
                                    gamma_scaling_factor,
                                    prevalence,
                                    mutvar_scaling_factor = 1,
                                    n_boot = 100) {
  mutvar_output_mod <- estimate_mutvar_trio(model,
                                              genetic_data,
                                              prevalence,
                                              gamma_scaling_factor = gamma_scaling_factor)



  bootstrap_mutvar_output_mod <- lapply(1:100,
                                              function(iter,model,genetic_data) {



                                                model_boot = model
                                                model_boot$conditional_likelihood = model_boot$conditional_likelihood[model$bootstrap_output$bootstrap_indices[,iter], , drop = FALSE]
                                                model_boot$features = model_boot$features[model$bootstrap_output$bootstrap_indices[,iter], , drop = FALSE]
                                                model_boot$delta = model$bootstrap_output$bootstrap_delta[[iter]]

                                                boot_mutvar = estimate_mutvar_trio(model = model_boot,
                                                                                               genetic_data = genetic_data[model$bootstrap_output$bootstrap_indices[,iter], , drop = FALSE],
                                                                                               prevalence = prevalence,
                                                                                               gamma_scaling_factor = gamma_scaling_factor)

                                              },
                                              model,
                                              genetic_data)

  bootstrap_mutvar_ests_mod = sapply(1:length(bootstrap_mutvar_output_mod), function(x) bootstrap_mutvar_output_mod[[x]]$total_mutvar)
  mutvar_CI_mod = quantile(bootstrap_mutvar_ests_mod,c(0.025,0.975))

  return(list(mutvar = mutvar_output_mod$total_mutvar * mutvar_scaling_factor,
              bootstrap_mutvar_ests = bootstrap_mutvar_ests_mod* mutvar_scaling_factor,
              mutvar_lower = mutvar_CI_mod[1]* mutvar_scaling_factor,
              mutvar_upper = mutvar_CI_mod[2]* mutvar_scaling_factor))

}


scaled_peneff <- function(model, genetic_data, prevalence, gamma_scaling_factor = 1) {
  peneff = effective_penetrance_func(model,
                                     genetic_data,
                                     prevalence,
                                     gamma_scaling_factor = gamma_scaling_factor)
  peneff_CI = NA

  cat("...bootstrap effective penetrance")
  bootstrap_peneff_ests = pbsapply(1:length(model$bootstrap_output$bootstrap_delta),
                                   function(iter) {
                                     model_boot = model
                                     model_boot$conditional_likelihood = model_boot$conditional_likelihood[model$bootstrap_output$bootstrap_indices[,iter], , drop = FALSE]
                                     model_boot$features = model_boot$features[model$bootstrap_output$bootstrap_indices[,iter], , drop = FALSE]
                                     model_boot$delta = model$bootstrap_output$bootstrap_delta[[iter]]

                                     genetic_data_boot = genetic_data[model$bootstrap_output$bootstrap_indices[,iter], , drop = FALSE]

                                     effective_penetrance_func(model_boot,
                                                               genetic_data_boot,
                                                               prevalence,
                                                               gamma_scaling_factor = gamma_scaling_factor)
                                   })

  peneff_CI = quantile(bootstrap_peneff_ests,c(0.025,0.975))

  peneff_df = data.frame(peneff = peneff,
                         peneff_lower = peneff_CI[1],
                         peneff_upper = peneff_CI[2],
                         RR_eff = peneff/prevalence,
                         RR_eff_lower = peneff_CI[1]/prevalence,
                         RR_eff_upper = peneff_CI[2]/prevalence)

  return(peneff_df)
}

bootstrap_mutvar <- function(model,genetic_data,prevalence, n_boot = 100) {
  bootstrap_mutvar_output <- pblapply(1:n_boot,
                                    function(iter) {



                                      model_boot = model
                                      model_boot$conditional_likelihood = model_boot$conditional_likelihood[model$bootstrap_output$bootstrap_indices[,iter], , drop = FALSE]
                                      model_boot$features = model_boot$features[model$bootstrap_output$bootstrap_indices[,iter], , drop = FALSE]
                                      model_boot$delta = model$bootstrap_output$bootstrap_delta[[iter]]

                                      boot_mutvar = estimate_mutvar_trio(model = model_boot,
                                                                         genetic_data = genetic_data[model$bootstrap_output$bootstrap_indices[,iter], , drop = FALSE],
                                                                         prevalence = prevalence)

                                    })

  bootstrap_mutvar_ests = sapply(1:length(bootstrap_mutvar_output), function(x) bootstrap_mutvar_output[[x]]$total_mutvar)
  mutvar_CI = quantile(bootstrap_mutvar_ests,c(0.025,0.975))

  return(list(ests = bootstrap_mutvar_ests,
              CI = mutvar_CI))
}

get_summarytable <- function(stratum_list,model_list,data,prev_factor) {
  summarytable <- data.frame()

  for (stratum in stratum_list) {
    PTV_index = which(data$loop_vars$names == paste0(stratum, " PTV"))
    Mis2_index = which(data$loop_vars$names == paste0(stratum, " Mis2"))

    mutvar_output_ptv = get_scaled_mutvar(model_list[[PTV_index]],
                                            gamma_scaling_factor = 1,
                                            genetic_data = get_genetic_data(PTV_index,data)$genetic_data,
                                            prevalence = data$loop_vars$prevalences[PTV_index] * prev_factor,
                                            mutvar_scaling_factor = data$ptv_scale_factor)

    mutvar_output_mis2 = get_scaled_mutvar(model_list[[Mis2_index]],
                                             gamma_scaling_factor = 1,
                                             genetic_data = get_genetic_data(Mis2_index,data)$genetic_data,
                                             prevalence = data$loop_vars$prevalences[Mis2_index] * prev_factor)

    mutvar_combined = mutvar_output_ptv$mutvar + mutvar_output_mis2$mutvar

    bootstrap_mutvar_combined = sapply(1:length(mutvar_output_ptv$bootstrap_mutvar_ests),
                                   function(x) mutvar_output_ptv$bootstrap_mutvar_ests[x] + mutvar_output_mis2$bootstrap_mutvar_ests[x])

    mutvar_combined_CI = quantile(bootstrap_mutvar_combined,c(0.025,0.975))


    combined_fraccase <- get_fraccase_df(data,model_list[[PTV_index]],
                                         get_genetic_data(PTV_index,data)$genetic_data,
                                         model_list[[Mis2_index]],
                                         get_genetic_data(Mis2_index,data)$genetic_data,
                                         gamma_scaling_factor = 1,
                                         RR_range = c(5))

    iter_df = data.frame(stratum = stratum,
                         prev = data$loop_vars$prevalences[PTV_index] * prev_factor,
                         mutvar_combined = mutvar_combined,
                         mutvar_combined_lower = mutvar_combined_CI[1],
                         mutvar_combined_upper = mutvar_combined_CI[2],
                         fraccase_RR5 = combined_fraccase$fraccases_greater_PTVplusMis2[combined_fraccase$RR_thresh == 5],
                         fraccase_RR5_lower = combined_fraccase$fraccases_greater_PTVplusMis2_lowerCI[combined_fraccase$RR_thresh == 5],
                         fraccase_RR5_upper = combined_fraccase$fraccases_greater_PTVplusMis2_upperCI[combined_fraccase$RR_thresh == 5],
                         PTV_peneff = model_list[[PTV_index]]$penetrance$effective_penetrance,
                         PTV_peneff_lower = model_list[[PTV_index]]$penetrance$effective_penetrance_CI[1],
                         PTV_peneff_upper = model_list[[PTV_index]]$penetrance$effective_penetrance_CI[2],
                         PTV_effRR = model_list[[PTV_index]]$penetrance$effective_penetrance/(data$loop_vars$prevalences[PTV_index] * prev_factor),
                         PTV_effRR_lower = model_list[[PTV_index]]$penetrance$effective_penetrance_CI[1]/(data$loop_vars$prevalences[PTV_index] * prev_factor),
                         PTV_effRR_upper = model_list[[PTV_index]]$penetrance$effective_penetrance_CI[2]/(data$loop_vars$prevalences[PTV_index] * prev_factor))

    summarytable = rbind(summarytable, iter_df)

  }
  return(summarytable)
}
