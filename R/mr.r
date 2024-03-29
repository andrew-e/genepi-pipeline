#'
#'
#'
#'
perform_mr_on_metabrain_datasets <- function(gwas_filename, ancestry="EUR", subcategory=NULL, exposures=c(), results_output) {
  file_pattern <- paste0(tolower(subcategory), "_", tolower(ancestry))
  metabrain_top_hits <- list.files(metabrain_top_hits_dir, pattern = file_pattern, full.names = T)

  run_mr_on_qtl_data(gwas_filename, results_output = results_output, qtl_files = metabrain_top_hits, exposures)
}


#' TODO: not in full use yet
#'
#'
perform_mr_on_pqtl_datasets <- function(gwas_filename, results_output, subcategory=None, cis_only = T) {
  pqtl_file_pattern <- if(cis_only) ".*_cis_only.*" else ".*_cis_trans.*"
  pqtl_datasets <- list.files(pqtl_top_hits_dir, pattern = pqtl_file_pattern, full.names = T)

  run_mr_on_qtl_data(gwas_filename,qtl_files = pqtl_datasets, results_output = results_output)
}

run_mr_on_qtl_data <- function(gwas_filename, qtl_files, results_output, exposures=c()) {
  all_qtl_mr_results <- lapply(qtl_files, function(qtl_file) {
    qtl_exposure <- TwoSampleMR::read_exposure_data(qtl_file,
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
    #qtl_exposure <- TwoSampleMR::clump_data(qtl_exposure, pop=ancestry)

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
    gwas_outcome_data$outcome <- file_prefix(qtl_file)

    harmonised_data <- TwoSampleMR::harmonise_data(qtl_exposure, gwas_outcome_data)
    mr_results <- TwoSampleMR::mr(harmonised_data, method_list=c("mr_wald_ratio", "mr_ivw"))
    if (length(exposures) > 0) {
      mr_results <- subset(mr_results, exposure %in% exposures)
    }

    mr_results$p.adjusted <- p.adjust(mr_results$pval, "fdr")

    qtl_dataset <- vroom::vroom(qtl_file) |>
      calculate_f_statistic()
    matching <- match(mr_results$exposure, qtl_dataset$EXPOSURE)
    mr_results$SNP <- qtl_dataset$SNP[matching]
    mr_results$F_STAT <- qtl_dataset$F_STAT[matching]

    return(mr_results)
  }) |> dplyr::bind_rows()


  vroom::vroom_write(all_qtl_mr_results, results_output)
}

#'
#'
compare_interesting_mr_results <- function(pqtl_mr_results, forest_plot_output_file) {
  pqtl_mr_results <- vroom::vroom(pqtl_mr_results)
  significant_pqtl_mr_results <- subset(pqtl_mr_results, p.adjusted < 0.05)

  interesting_pqtl_mr_results <- subset(pqtl_mr_results, exposure %in% significant_pqtl_mr_results$exposure) %>%
    dplyr::select(exposure, dplyr::everything())
  interesting_pqtl_mr_results$BETA <- interesting_pqtl_mr_results$b
  interesting_pqtl_mr_results$SE <- interesting_pqtl_mr_results$se

  grouped_forest_plot(interesting_pqtl_mr_results,
                      "Comparing MR Results across pQTL data sets",
                      "outcome",
                      forest_plot_output_file
  )
}
