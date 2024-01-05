
populate_gene_ids <- function(gwas) {
  if ("ENSEMBL_ID" %in% colnames(gwas) && !"GENE_ID" %in% colnames(gwas)) {
    return(ensembl_to_gene())
  }
  else if ("GENE_ID" %in% colnames(gwas) && !"ENSEMBL_ID" %in% colnames(gwas)) {
    return(gene_to_ensembl())
  }
  return(gwas)
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

populate_rsid_options <- list(no="NO", partial="PARTIAL", full="FULL")
rsid_builds <- list(GRCh37="b37_dbsnp156")

populate_rsid <- function(gwas, option=populate_rsid_options$partial) {
  if (!option %in% populate_rsid_options) stop(paste("Error: invalid option:", option))

  start <- Sys.time()
  if (option == populate_rsid_options$no || column_map$default$RSID %in% colnames(gwas)) {
    message("Skipping RSID population for GWAS")
  }
  else if (option == populate_rsid_options$full) {
    gwas <- populate_full_rsids(gwas)
  }
  else if (option == populate_rsid_options$partial) {
    gwas <- populate_partial_rsids(gwas)
  }
  message("RSID population time taken:")
  print(Sys.time() - start)
  return (gwas)
}

populate_partial_rsids <- function(gwas) {
  message("populating RSIDs based on 1000genomes...")
  marker_to_rsid_file <- paste0(genomic_data_dir, "1000genomes/marker_to_rsid.tsv.gz")
  chrpos_to_rsid <- vroom::vroom(marker_to_rsid_file, col_select=c("HG37", "RSID"))
  gwas$RSID <- chrpos_to_rsid$RSID[match(gwas$SNP, chrpos_to_rsid$HG37)]

  return(gwas)
}

populate_full_rsids <- function(gwas, build=rsid_builds$GRCh37) {
  dbsnp_dir <- paste0(genomic_data_dir, "dbsnp")
  if (!build %in% rsid_builds) stop(paste("Error: invalid rsid build option:", build))
  gc()

  gwas <- data.table::as.data.table(gwas)
  future::plan(future::multisession, workers = number_of_cpus_available)
  gwas <- genepi.utils::chrpos_to_rsid(gwas, "CHR", "BP", "EA", "OA", flip = "allow", dbsnp_dir=dbsnp_dir, build=build, alt_rsids = F)
  gwas <- tibble::as_tibble(gwas)
  return(gwas)
}
