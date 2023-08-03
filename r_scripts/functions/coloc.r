
#' run_coloc_analysis takes two already harmonised gwases, and runs coloc on the results
#' @param first_gwas: first gwas to be run through coloc.  This is the gwas that results will be based off
#' @param second_gwas: second gwas to be run through coloc
#' @returns
coloc_analysis <- function(first_gwas, second_gwas, output_file) {
    library(vroom, quietly=T)
    library(coloc, quietly=T)
    create_dir_for_files(output_file)

    coloc_columns <- c("SNP", "P", "SE", "N", "EAF")

    if (!is.data.frame(first_gwas)) {
        first_gwas <- vroom::vroom(first_gwas, col_select = coloc_columns)
    }
    if (!is.data.frame(second_gwas)) {
        second_gwas <- vroom::vroom(second_gwas, col_select = coloc_columns)
    }

    if (!("N" %in% colnames(first_gwas$N)) | !("N" %in% colnames(second_gwas$N))) {
        stop("GWASes must contain column N (sample size) to complete coloc analysis")
    }

    first_gwas$P[first_gwas$P == 0] <- .Machine$double.xmin
    second_gwas$P[second_gwas$P == 0] <- .Machine$double.xmin

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
      ~exposure, ~coloc_h0, ~coloc_h1, ~coloc_h2, ~coloc_h3, ~coloc_h4,
      "huh", result$summary[2], result$summary[3], result$summary[4], result$summary[5], result$summary[6]
    )

    vroom::vroom_write(coloc_results, output_file)

    return(coloc_results)
}
