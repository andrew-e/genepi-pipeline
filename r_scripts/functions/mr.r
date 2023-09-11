#'
#'
perform_mr_on_pqtl_datasets <- function(gwas_filename, results_output, cis_only = T, pop="EUR") {
  library(dplyr)
  library(vroom)
  library(TwoSampleMR)

  pqtl_file_pattern <- if(cis_only) ".*_cis_only.*" else ".*_cis_trans.*"
  pqtl_datasets <- list.files(pqtl_top_hits_dir, pattern = pqtl_file_pattern, full.names = T)

  all_pqtl_mr_results <- lapply(pqtl_datasets, function(pqtl_dataset) {
     pqtl_exposure <- TwoSampleMR::read_exposure_data(pqtl_dataset,
                                                      snp_col="SNP",
                                                      effect_allele_col="EA",
                                                      other_allele_col="OA",
                                                      eaf_col="EAF",
                                                      beta_col="BETA",
                                                      se_col="SE",
                                                      pval_col="P",
                                                      phenotype_col="EXPOSURE",
                                                      sep="\t"
    )

    #Don't need to clump only because the pqtl cis only data sets are already clummped by ld
    #pqtl_exposure <- TwoSampleMR::clump_data(pqtl_exposure, pop=pop)

    gwas_outcome_data <- TwoSampleMR::read_outcome_data(gwas_filename,
                                                        snp_col = "SNP",
                                                        beta_col = "BETA",
                                                        se_col = "SE",
                                                        effect_allele_col = "EA",
                                                        other_allele_col = "OA",
                                                        pval_col = "P",
                                                        eaf_col = "EAF",
                                                        snps = pqtl_exposure$SNP,
                                                        sep="\t",
    )
    gwas_outcome_data$outcome <- file_prefix(pqtl_dataset)

    harmonised_data <- TwoSampleMR::harmonise_data(pqtl_exposure, gwas_outcome_data)
    mr_results <- TwoSampleMR::mr(harmonised_data, method_list=c("mr_wald_ratio", "mr_ivw"))
    mr_results$p.adjusted <- p.adjust(mr_results$pval, "fdr")

    return(mr_results)
  }) %>% dplyr::bind_rows()

  vroom::vroom_write(all_pqtl_mr_results, results_output)
}

#'
#'
compare_interesting_mr_results <- function(pqtl_mr_results, forest_plot_output_file) {
  library(vroom)
  library(dplyr)
  pqtl_mr_results <- vroom::vroom(pqtl_mr_results)
  significant_pqtl_mr_results <- subset(pqtl_mr_results, p.adjusted < 0.05)

  interesting_pqtl_mr_results <- subset(pqtl_mr_results, exposure %in% significant_pqtl_mr_results$exposure) %>%
    dplyr::select(exposure, everything())
  interesting_pqtl_mr_results$BETA <- interesting_pqtl_mr_results$b
  interesting_pqtl_mr_results$SE <- interesting_pqtl_mr_results$se

  grouped_forest_plot(interesting_pqtl_mr_results,
                      "Comparing MR Results across pQTL data sets",
                      "outcome",
                      forest_plot_output_file
  )
}