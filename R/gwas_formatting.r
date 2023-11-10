library(vroom, quietly=T)
library(tidyr, quietly=T)
library(dplyr, quietly=T)

column_map <- list(
  default = list(SNP="SNP", CHR="CHR", BP="BP", EA="EA", OA="OA", EAF="EAF", P="P", Z="Z", BETA="BETA", SE="SE", OR="OR", OR_SE="OR_SE", OR_LB="OR_LB", OR_UB="OR_UB", RSID="RSID", N="N"),
  metal = list(SNP="MarkerName", EA="Allele1", OA="Allele2", EAF="Freq1", P="P-value", BETA="Effect", SE="StdErr"),
  ieu_ukb = list(SNP="SNP", BETA="BETA", SE="SE", EA="ALLELE1", OA="ALLELE0", EAF="A1FREQ", P="P_BOLT_LMM_INF")
)

standardise_gwas <- function(file_gwas,
                             output_file,
                             input_format="default",
                             populate_rsid_option="PARTIAL",
                             bespoke_column_map=NULL) {

  if (is.null(column_map[[input_format]])) {
    stop(paste("Error: invalid input_format!", input_format, "is not recognised."))
  }

  #TODO: if we need to add bespoke input format wrangling here, we can
  gwas_columns <- if(!is.null(bespoke_column_map)) bespoke_column_map else column_map[[input_format]]

  gwas <- vroom::vroom(file_gwas) |>
    change_column_names(gwas_columns) |>
    filter_incomplete_rows() |>
    standardise_columns() |>
    standardise_alleles() |>
    health_check() |>
    populate_rsid(populate_rsid_option)

  vroom::vroom_write(gwas, output_file)
}

populate_rsid <- function(option) {
  if (option == "NO" || column_map$default$RSID %in% colnames(gwas)) {
    message("Skipping RSID population for GWAS")
  }
  else if (option == "FULL") {
    gwas <- populate_full_rsids(gwas)
  }
  else if (option == "PARTIAL") {
    gwas <- populate_partial_rsids(gwas)
  }
  return (gwas)
}

format_gwas_output <- function(file_gwas, output_file, output_format="default") {
  gwas <- vroom::vroom(file_gwas) |>
    change_column_names(column_map[[output_format]], opposite_mapping = T)

  vroom::vroom_write(gwas, output_file)
}

filter_incomplete_rows <- function(gwas) {
  filtered_gwas <- gwas[!is.na(gwas$EAF) & !is.null(gwas$EAF) &
                          !is.na(gwas$OA) & !is.null(gwas$OA) &
                          !is.na(gwas$EA) & !is.null(gwas$EA) &
                          !is.na(gwas$CHR) & !is.null(gwas$CHR) &
                          !is.na(gwas$BP) & !is.null(gwas$BP),
  ]
  filtered_rows <- nrow(gwas) - nrow(filtered_gwas)
  if (filtered_rows > 0) {
    warning(paste("Warning: Filtering out ", filtered_rows, "rows due to NULLs and NAs"))
  }
  return(filtered_gwas)
}

standardise_columns <- function(gwas) {
  if (!all(c("CHR", "BP") %in% colnames(gwas))) {
    if(all(grepl("\\d:\\d", gwas$SNP))) {
      gwas <- tidyr::separate(data = gwas, col = "SNP", into = c("CHR", "BP"), sep = "[:_]", remove = F)
    }
  }

  if (all(c("OR", "OR_LB", "OR_UB") %in% colnames(gwas)) & !all(c("BETA", "SE") %in% colnames(gwas))) {
    gwas <- convert_or_to_beta(gwas)
  }
  else if (all(c("OR", "OR_SE") %in% colnames(gwas)) & !all(c("BETA", "SE") %in% colnames(gwas))) {
    gwas <- convert_or_to_beta(gwas)
  }

  gwas$BP <- as.numeric(gwas$BP)
  gwas$P <- as.numeric(gwas$P)
  gwas$P[gwas$P == 0] <- .Machine$double.xmin

  return(gwas)
}



health_check <- function(gwas) {
  if (nrow(gwas[gwas$P <= 0 | gwas$P > 1, ]) > 0) {
    stop("GWAS has some P values outside accepted range")
  }
  if (nrow(gwas[gwas$EAF < 0 | gwas$EAF > 1, ]) > 0) {
    stop("GWAS has some EAF values outside accepted range")
  }
  #if ("OR" %in% colnames(gwas) & nrow(gwas[gwas$OR < 0, ]) > 0) {
  #  stop("GWAS has some OR values outside accepted range")
  #}
  return(gwas)
}


#' change_column_names: private function that takes a named list of column names
#'  and changes the supplied data frame's column names accordingly
#'
#' @param: gwas dataframe of gwas to standardise column names
#' @param: columns named list for
#' @param: opposite_mapping logical flag on if we are mapping from key to value or vice verca
change_column_names <- function(gwas, columns = list(), opposite_mapping = FALSE) {
  if (!opposite_mapping) {
    for (name in names(columns)) {
      names(gwas)[names(gwas) == columns[name]] <- name
    }
  }
  else {
    for (name in names(columns)) {
      names(gwas)[names(gwas) == name] <- columns[name]
    }
  }
  return(gwas)
}

standardise_alleles <- function(gwas) {
  gwas$EA <- toupper(gwas$EA)
  gwas$OA <- toupper(gwas$OA)

  to_flip <- gwas$EA > gwas$OA
  gwas$EAF[to_flip] <- 1 - gwas$EAF[to_flip]
  gwas$BETA[to_flip] <- -1 * gwas$BETA[to_flip]

  temp <- gwas$OA[to_flip]
  gwas$OA[to_flip] <- gwas$EA[to_flip]
  gwas$EA[to_flip] <- temp

  gwas$SNP <- toupper(paste0(gwas$CHR, ":", format(gwas$BP, scientific = F, trim = T), "_", gwas$EA, "_", gwas$OA))
  gwas <- dplyr::select(gwas, SNP, CHR, BP, EA, OA, EAF, BETA, SE, P, dplyr::everything())

  return(gwas)
}

#' harmonise_gwases: takes a list of gwases, get the SNPs in common
#' across all datasets arranged to be in the same order
#'
#' @param: elipses of gwases
#' @return: list of harmonised gwases
#'
harmonise_gwases <- function(...) {
  gwases <- list(...)

  snpids <- Reduce(intersect, lapply(gwases, function(gwas) gwas$SNP))
  message(paste("Number of shared SNPs after harmonisation:", length(snpids)))

  gwases <- lapply(gwases, function(gwas) {
    dplyr::filter(gwas, SNP %in% snpids) |>
    dplyr::arrange(SNP)
  })

  return(gwases)
}

populate_snp_from_rsid <- function(gwas) {
  start <- Sys.time()
  marker_to_rsid_file <- paste0(thousand_genomes_dir, "marker_to_rsid_full.tsv.gz")
  marker_to_rsid <- vroom::vroom(marker_to_rsid_file, col_select=c("HG37", "RSID"))
  message(paste("loaded file: ", Sys.time()-start))

  matching <- match(gwas$RSID, marker_to_rsid$RSID)
  gwas$CHRBP <- marker_to_rsid$HG37[matching]
  gwas <- tidyr::separate(data = gwas, col = "CHRBP", into = c("CHR", "BP"), sep = ":")
  message(paste("mapped and returned: ", Sys.time()-start))

  return(gwas)
}



create_bed_file_from_gwas <- function(gwas, output_file) {
  N <- nrow(gwas)
  bed_file <- data.frame(CHR = character(N),
    BP1 = numeric(N),
    BP2 = numeric(N)
  )

  split <- tidyr::separate(data = gwas, col = "SNP", into = c("CHR", "BP1"), sep = "[:_]", remove = T)

  bed_file$CHR <- paste0("chr", split$CHR)
  bed_file$BP1 <- split$BP1
  bed_file$BP2 <- split$BP1

  vroom::vroom_write(bed_file, output_file, col_names=F)
}

read_liftover_output_file <- function(filename) {
  liftover_bed <- vroom::vroom(filename)
  N <- nrow(liftover_bed)

  bed_map <- data.frame(
    ORIG_MARKER = character(N),
    NEW_MARKER = character(N)
  )

  liftover_bed$V1 <- gsub("chr", "", liftover_bed$V1)
  liftover_bed$V4 <- gsub("^.*-", "", liftover_bed$V4)
  bed_map$NEW_MARKER <- paste(liftover_bed$V1, liftover_bed$V2, sep=":")
  bed_map$ORIG_MARKER <- paste(liftover_bed$V1, liftover_bed$V4, sep=":")

  return(bed_map)
}