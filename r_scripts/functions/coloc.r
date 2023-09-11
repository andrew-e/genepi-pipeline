#' run_coloc_analysis takes two already harmonised gwases, and runs coloc on the results
#' @param first_gwas: first gwas to be run through coloc.  This is the gwas that results will be based off
#' @param second_gwas: second gwas to be run through coloc
#' @returns
coloc_analysis <- function(first_gwas_filename, second_gwas_filename, output_file, exposure_name, chr=NA, bp=NA, range=NA) {
  library(dplyr, quietly=T)
  library(vroom, quietly=T)
  library(coloc, quietly=T)

  coloc_columns <- c("SNP", "CHR", "BP", "P", "SE", "N", "EAF")
  #if (!is.na(chr)) {
  #  a=proc.time(); first_gwas <- vroom_chr(first_gwas_filename, chr, coloc_columns); b=proc.time(); b-a
  #  second_gwas <- vroom_chr(second_gwas_filename, chr, coloc_columns)
  #}
  #else {
    first_gwas <- vroom::vroom(first_gwas_filename, col_select = coloc_columns)
    second_gwas <- vroom::vroom(second_gwas_filename, col_select = coloc_columns)
  #}

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

  #vroom::vroom_write(coloc_results, output_file)
  return(coloc_results)
}
