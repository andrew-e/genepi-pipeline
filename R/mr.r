perform_mr_on_pqtl_datasets <- function(gwas_filename, ancestry="EUR", results_output, cis_only = T) {
  pqtl_file_pattern <- if(cis_only) ".*_cis_only.*" else ".*_cis_trans.*"
  pqtl_datasets <- list.files(pqtl_top_hits_dir, pattern = pqtl_file_pattern, full.names = T)

  generic_mr_qtl_comparison(gwas_filename, results_output = results_output, qtl_datasets = pqtl_datasets, pop = ancestry)
}

perform_mr_on_metabrain_datasets <- function(gwas_filename, ancestry="EUR", results_output) {
  file_pattern <- paste0("_", tolower(ancestry))
  qtl_datasets <- list.files(metabrain_top_hits_dir, pattern = file_pattern, full.names = T)

  generic_mr_qtl_comparison(gwas_filename, results_output = results_output, qtl_datasets = qtl_datasets, pop = ancestry)
}


generic_mr_qtl_comparison <- function(gwas_filename, qtl_datasets, results_output, ancestry="EUR") {
  all_pqtl_mr_results <- lapply(qtl_datasets, function(qtl_dataset) {
    qtl_exposure <- TwoSampleMR::read_exposure_data(qtl_dataset,
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

    #Don't need to clump because the qtl data sets are already only top hits
    #qtl_exposure <- TwoSampleMR::clump_data(pqtl_exposure, pop=ancestry)

    gwas_outcome_data <- TwoSampleMR::read_outcome_data(gwas_filename,
                                                        snp_col = "SNP",
                                                        beta_col = "BETA",
                                                        se_col = "SE",
                                                        effect_allele_col = "EA",
                                                        other_allele_col = "OA",
                                                        pval_col = "P",
                                                        eaf_col = "EAF",
                                                        snps = qtl_exposure$SNP,
                                                        sep="\t",
    )
    gwas_outcome_data$outcome <- file_prefix(qtl_dataset)

    harmonised_data <- TwoSampleMR::harmonise_data(qtl_exposure, gwas_outcome_data)
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
