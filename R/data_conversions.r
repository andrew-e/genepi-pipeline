#' convert_or_to_beta: Given an OR and lower and upper bounds,
#'   calculates the BETA, and SE.
#'   based on this answer: https://stats.stackexchange.com/a/327684
#'
#' @param gwas: dataframe with the following columns: OR, LB (lower bound), UB (upper bound)
#' @return gwas with new columns BETA and SE
#'
convert_or_to_beta <- function(gwas) {
  gwas$BETA <- log(gwas$OR)

  if ("OR_SE" %in% colnames(gwas)) {
    gwas$OR_UB <- gwas$OR + (gwas$OR_SE * 1.96)
    gwas$OR_LB <- gwas$OR - (gwas$OR_SE * 1.96)
    gwas$SE <- gwas$OR_SE
  }
  if (all(c("OR_LB", "OR_UB") %in% colnames(gwas))) {
    gwas$SE <- (gwas$OR_UB - gwas$OR_LB) / (2 * 1.96)
  }
  else {
    stop("Need OR_SE, or OR_LB + OR_UB to complete conversion")
  }

  return(gwas)
}

convert_beta_to_or <- function(gwas) {
  gwas$OR <- exp(gwas$BETA)
  gwas$OR_SE <- gwas$SE
  return(gwas)
}

convert_z_to_p <- function(gwas) {
  gwas$P <- 2 * pnorm(-abs(gwas$Z))
  return(gwas)
}

calculate_f_statistic <- function(gwas) {
  gwas$F_STAT <- qchisq(gwas$P,1,low=F)
  return(gwas)
}

populate_gene_ids <- function(gwas) {
  if ("ENSEMBL_ID" %in% colnames(gwas) && !"GENE_ID" %in% colnames(gwas)) {
    ensembl_to_gene()
  }
  else if ("GENE_ID" %in% colnames(gwas) && !"ENSEMBL_ID" %in% colnames(gwas)) {
    gene_to_ensembl()
  }
}

ensembl_to_gene <- function(gwas) {
  library(EnsDb.Hsapiens.v79, lib.loc=paste0(genomic_data_dir, "ensembl/"))
  gwas <- ensembldb::select(EnsDb.Hsapiens.v79, keys = gwas$ENSEMBL_ID, keytype = "GENEID", columns = c("SYMBOL","GENEID")) |>
      dplyr::rename(GENE_ID = "GENEID", ENSEMBL_ID = "SYMBOL") |>
      dplyr::merge(gwas, by="ENSEMBL_ID")
  return(gwas)
}

gene_to_ensembl <- function(gwas) {
  library(EnsDb.Hsapiens.v79, lib.loc=paste0(genomic_data_dir, "ensembl/"))
  gwas <- ensembldb::select(EnsDb.Hsapiens.v79, keys = gwas$GENE_ID, keytype = "SYMBOL", columns = c("SYMBOL","GENEID")) |>
      dplyr::rename(GENE_ID = "GENEID", ENSEMBL_ID = "SYMBOL") |>
      dplyr::merge(gwas, by="ENSEMBL_ID")
  return(gwas)
}

populate_rsid <- function(gwas, option) {
  message("Start Time:")
  message(Sys.time())
  if (option == "NO" || column_map$default$RSID %in% colnames(gwas)) {
    message("Skipping RSID population for GWAS")
  }
  else if (option == "FULL") {
    gwas <- populate_full_rsids(gwas)
  }
  else if (option == "PARTIAL") {
    gwas <- populate_partial_rsids(gwas)
  }
  message("Time taken:")
  message(Sys.time())
  return (gwas)
}

populate_partial_rsids <- function(gwas) {
  message("populating RSIDs based on 1000genomes...")
  marker_to_rsid_file <- paste0(genomic_data_dir, "1000genomes/marker_to_rsid.tsv.gz")
  chrpos_to_rsid <- vroom::vroom(marker_to_rsid_file, col_select=c("HG37", "RSID"))
  gwas$RSID <- chrpos_to_rsid$RSID[match(gwas$SNP, chrpos_to_rsid$HG37)]

  return(gwas)
}

populate_full_rsids <- function(gwas, build="b37_dbsnp156") {
  gc()
  gwas <- data.table::as.data.table(gwas)
  dbsnp_dir <- paste0(genomic_data_dir, "dbsnp")
  future::plan(future::multisession, workers = number_of_cpus_available)
  gwas <- genepi.utils::chrpos_to_rsid(gwas, "CHR", "BP", "EA", "OA", flip = "allow", dbsnp_dir=dbsnp_dir, build=build, alt_rsids = F)
  gwas <- tibble::as_tibble(gwas)
  return(gwas)
}