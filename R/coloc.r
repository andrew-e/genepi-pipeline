run_coloc_on_qtl_mr_results <- function(mr_results_file,
                                        gwas_file,
                                        qtl_dataset,
                                        exposures=c(),
                                        output_file) {
  coloc_columns <- c("SNP", "CHR", "BP", "P", "SE", "N", "EAF")
  gwas <- vroom::vroom(gwas_file, col_select = coloc_columns) |> subset(EAF > 0 & EAF < 1)
  range <- 500000

  mr_results <- vroom::vroom(mr_results_file)
  if (length(exposures) > 0) {
    mr_results <- subset(mr_results, exposure %in% exposures)
  }
  else {
    mr_results <- head(mr_results[order(mr_results$p.adjusted), ], 10)
  }
  mr_results <- tidyr::separate(data = mr_results, col = "SNP", into = c("CHR", "BP"), sep = "[:_]", remove = F)

  coloc_results <- apply(mr_results, 1, function(mr_result) {
    if (qtl_dataset == "metabrain") {
      outcome <- unlist(strsplit(mr_result[["outcome"]], "_"))
      brain_region <- outcome[1]
      ancestry <- outcome[2]

      chr <- as.numeric(mr_result[["CHR"]])
      bp <- as.numeric(mr_result[["BP"]])

      metabrain_dir <- "scratch/data/qtls/metabrain/"
      qtl_gwas_file <- paste0(metabrain_dir, "gwas/", brain_region, "/", mr_result[["exposure"]], "_", ancestry, ".tsv.gz")

      qtl_chr_gwas <- vroom::vroom(qtl_gwas_file) |> subset(EAF > 0 & EAF < 1)
    }
    else {
      stop("Error: qtl dataset not supported for coloc right now")
    }

    result <- coloc_analysis(gwas, qtl_chr_gwas, mr_result[["exposure"]], chr, bp, range)

    miami_filename <- paste0("scratch/results/plots/mr_metabrain_coloc_", file_prefix(gwas_file), "_", mr_result[["exposure"]], ".png")
    miami_plot(qtl_gwas_file, gwas_file, miami_filename, paste0("Miami Plot of", mr_result[["exposure"]]), chr, bp, range)

    return(result)
  }) |> dplyr::bind_rows()

  vroom::vroom_write(coloc_results, output_file)
}

#' run_coloc_analysis takes two already harmonised gwases, and runs coloc on the results
#' @param first_gwas: first gwas to be run through coloc.  This is the gwas that results will be based off
#' @param second_gwas: second gwas to be run through coloc
#' @returns
coloc_analysis <- function(first_gwas, second_gwas, exposure_name, chr=NA, bp=NA, range=NA) {
  library(dplyr, quietly=T)
  library(coloc, quietly=T)

  if (!is.na(chr) & !is.na(bp) & !is.na(range)) {
    first_gwas <- gwas_region(first_gwas, chr, bp, range)
    second_gwas <- gwas_region(second_gwas, chr, bp, range)
  }

  harmonised_gwases <- harmonise_gwases(first_gwas, second_gwas)
  first_gwas <- harmonised_gwases[[1]]
  second_gwas <- harmonised_gwases[[2]]

  if (!("N" %in% colnames(first_gwas)) | !("N" %in% colnames(second_gwas))) {
    stop("GWASes must contain column N (sample size) to complete coloc analysis")
  }

  first_coloc_dataset <- list(
    pvalues = first_gwas$P,
    N = first_gwas$N[1],
    varbeta = first_gwas$SE^2,
    type = "quant",
    snp = first_gwas$SNP,
    MAF = first_gwas$EAF
  )

  second_coloc_dataset <- list(
    pvalues = second_gwas$P,
    N = second_gwas$N[1],
    varbeta = second_gwas$SE^2,
    type = "quant",
    snp = second_gwas$SNP,
    MAF = second_gwas$EAF
  )

  result <- coloc::coloc.abf(dataset1 = first_coloc_dataset, dataset2 = second_coloc_dataset)
  coloc_results <- tibble::tribble(
    ~exposure, ~h0, ~h1, ~h2, ~h3, ~h4,
    exposure_name, result$summary[2], result$summary[3], result$summary[4], result$summary[5], result$summary[6]
  )

  return(coloc_results)
}
