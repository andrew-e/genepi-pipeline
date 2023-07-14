coloc_analysis <- function(first_gwas, second_gwas, chr, bp, range = 250000) {
    suppressPackageStartupMessages(library(coloc))
    range = floor(range / 2)

    if (!is.data.frame(first_gwas)) {
        first_gwas <- data.table::fread(first_gwas)
    }
    if (!is.data.frame(second_gwas)) {
        second_gwas <- data.table::fread(second_gwas)
    }

    first_gwas$P[first_gwas$P == 0] <- .Machine$double.xmin
    second_gwas$P[second_gwas$P == 0] <- .Machine$double.xmin

    harmonised <- harmonise_gwases(first_gwas, second_gwas, chr=chr, bp=bp, range=range)
    first_gwas <- harmonised[[1]]
    second_gwas <- harmonised[[2]]
    
    miami_plot(first_gwas, second_gwas)

    first_coloc_dataset <- list(
        pvalues = first_gwas$P,
        N = first_gwas$N[1],
        varbeta = first_gwas$SE^2,
        type = "quant",
        snp = first_gwas$MARKER,
        MAF = first_gwas$EAF
    )

    second_coloc_dataset = list(
        pvalues = second_gwas$P,
        N = second_gwas$N[1],
        varbeta = second_gwas$SE^2,
        type = "quant",
        snp = second_gwas$MARKER,
        MAF = second_gwas$EAF
    )
    
    result <- coloc.abf(dataset1 = first_coloc_dataset, dataset2 = second_coloc_dataset)
    return(result)
}
