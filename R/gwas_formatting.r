library(vroom, quietly=T)
library(tidyr, quietly=T)
library(dplyr, quietly=T)

column_map <- list(
  default = list(SNP="SNP", CHR="CHR", BP="BP", EA="EA", OA="OA", EAF="EAF", P="P", BETA="BETA", SE="SE", OR="OR", OR_SE="OR_SE", OR_LB="OR_LB", OR_UB="OR_UB", RSID="RSID", N="N"),
  metal = list(SNP="MarkerName", EA="Allele1", OA="Allele2", EAF="Freq1", P="P-value", BETA="Effect", SE="StdErr"),
  ieu_ukb = list(SNP="SNP", BETA="BETA", SE="SE", EA="ALLELE1", OA="ALLELE0", EAF="A1FREQ", P="P_BOLT_LMM_INF")
)

standardise_gwas <- function(file_gwas,
                             output_file,
                             input_format="default",
                             populate_rsid=F,
                             bespoke_column_map=NULL) {

  if (is.null(column_map[[input_format]])) {
    stop(paste("Error: invalid input_format!", input_format, "is not recognised."))
  }

  #TODO: if we need to add bespoke input format wrangling here, we can
  gwas_columns <- if(!is.null(bespoke_column_map)) bespoke_column_map else column_map[[input_format]]

  gwas <- vroom::vroom(file_gwas) %>%
    change_column_names(gwas_columns) %>%
    filter_incomplete_rows() %>%
    standardise_columns() %>%
    standardise_alleles() %>%
    health_check() %>%
    populate_rsid_from_1000_genomes(populate_rsid)

  vroom::vroom_write(gwas, output_file)
}

format_gwas_output <- function(file_gwas, output_file, output_format="default") {
  gwas <- vroom::vroom(file_gwas) %>%
    change_column_names(column_map[[output_format]], opposite_mapping = T)

  vroom:vroom_write(gwas, output_file, delim="\t")
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
    print(paste("WARNING: Filtering out ", nrow(gwas) - nrow(filtered_gwas), "rows due to NULLs and NAs"))
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


#' convert_or_to_beta: Given an OR and lower and upper bounds,
#'   calculates the BETA, and SE.
#'   based on this answer: https://stats.stackexchange.com/a/327684
#'
#' @param gwas: dataframe with the following columns: OR, LB (lower bound), UB (upper bound)
#' @return gwas with new columns BETA and SE
#'
convert_or_to_beta <- function(gwas) {
  library(boot)
  gwas$BETA <- log(gwas$OR)

  if ("OR_SE" %in% colnames(gwas)) {
    gwas$OR_UB <- gwas$OR + (gwas$OR_SE * 1.96)
    gwas$OR_LB <- gwas$OR - (gwas$OR_SE * 1.96)
  }

  gwas$SE <- (gwas$OR_UB - gwas$OR_LB) / (2 * 1.95996)
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
  gwas %>% dplyr::select(SNP, CHR, BP, EA, OA, EAF, BETA, SE, P, everything())

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
  print(paste("Number of shared SNPs after harmonisation:", length(snpids)))

  gwases <- lapply(gwases, function(gwas) {
    gwas %>%
      dplyr::filter(SNP %in% snpids) %>%
      dplyr::arrange(SNP)
  })

  return(gwases)
}

populate_snp_from_rsid <- function(gwas) {
  start <- Sys.time()
  marker_to_rsid_file <- paste0(thousand_genomes_dir, "marker_to_rsid_full.tsv.gz")
  marker_to_rsid <- vroom::vroom(marker_to_rsid_file, col_select=c("HG37", "RSID"))
  print(paste("loaded file: ", Sys.time()-start))

  matching <- match(gwas$RSID, marker_to_rsid$RSID)
  gwas$CHRBP <- marker_to_rsid$HG37[matching]
  gwas <- tidyr::separate(data = gwas, col = "CHRBP", into = c("CHR", "BP"), sep = ":")
  print(paste("mapped and returned: ", Sys.time()-start))

  return(gwas)
}

populate_rsid_from_1000_genomes <- function(gwas, populate_rsid=F) {
  if (populate_rsid == F) return(gwas)

  #if (populate_rsid == "full") {
  #  marker_to_rsid_file <- paste0(thousand_genomes_dir, "marker_to_rsid_full.tsv.gz")
  #}
  #else {
  #  marker_to_rsid_file <- paste0(thousand_genomes_dir, "marker_to_rsid.tsv.gz")
  #}

  if(column_map$default$RSID %in% colnames(gwas)) {
    print("GWAS already has an RSID field, will not overwrite")
    return(gwas)
  }
  print("populating RSID...")
  marker_to_rsid_file <- paste0(genome_data_dir, "marker_to_rsid.tsv.gz")
  chrpos_to_rsid <- vroom::vroom(marker_to_rsid_file, col_select=c("HG37", "RSID"))
  gwas$RSID <- chrpos_to_rsid$RSID[match(gwas$SNP, chrpos_to_rsid$HG37)]

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

  vroom::vroom_write(bed_file, output_file)
}