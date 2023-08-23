library(vroom, quietly=T)
library(tidyr, quietly=T)
library(dplyr, quietly=T)

column_map <- list(
  default = list(SNP="SNP", EA="EA", OA="OA", EAF="EAF", P="P", BETA="BETA", SE="SE", OR="OR", OR_LB="OR_LB", OR_UB="OR_UB", RSID="RSID"),
  metal = list(SNP="MarkerName", EA="Allele1", OA="Allele2", EAF="Freq1", P="P-value", BETA="Effect", SE="StdErr"),
  ieu_ukb = list(SNP="SNP", BETA="BETA", SE="SE", EA="ALLELE1", OA="ALLELE0", EAF="A1FREQ", P="P_BOLT_LMM_INF")
)

standardise_gwas <- function(file_gwas, output_file, input_format="default", populate_rsid=NULL) {
  if (is.null(column_map[[input_format]])) {
    stop(paste("Error: invalid input_format!", input_format, "is not recognised."))
  }
  create_dir_for_files(output_file)
  
  #TODO: if we need to add bespoke input format wrangling here, we can
  
  gwas <- vroom::vroom(file_gwas) %>%
    change_column_names(column_map[[input_format]]) %>%
    standardise_alleles() %>%
    standardise_columns() %>%
    health_check() %>%
    populate_rsid_from_1000_genomes(populate_rsid)

  vroom::vroom_write(gwas, output_file)
}

format_gwas_output <- function(file_gwas, output_file, output_format="default") {
  gwas <- vroom::vroom(file_gwas) %>%
    change_column_names(column_map[[output_format]], opposite_mapping = T)

  vroom:vroom_write(gwas, output_file, delim="\t")
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

  gwas$SNP <- toupper(gwas$SNP)
  gwas$SNP[grep("^RS", gwas$SNP)] <- tolower(gwas$SNP)
  gwas$EA <- toupper(gwas$EA)
  gwas$OA <- toupper(gwas$OA)

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

    gwas$SE <- (gwas$OR_UB - gwas$OR_LB) / (2 * 1.95996)
    return(gwas)
}

health_check <- function(gwas) {
  if (nrow(gwas[gwas$P < 0 | gwas$P > 1, ]) > 0) {
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

  gwas$SNP <- paste0(gwas$CHR, ":", format(gwas$BP, scientific = F, trim = T), "_", gwas$EA, "_", gwas$OA)
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

populate_rsid_from_1000_genomes <- function(gwas, populate_rsid=NULL, reference_build="hg37") {
  if (is.null(populate_rsid)) return(gwas)

  if (populate_rsid == "full") {
    marker_to_rsid_file <- paste0(thousand_genomes_dir, "marker_to_rsid_full.tsv.gz")
  }
  else {
    marker_to_rsid_file <- paste0(thousand_genomes_dir, "marker_to_rsid.tsv.gz")
  }

  if(column_map$default$RSID %in% colnames(gwas)) {
    print("GWAS already has an RSID field, will not overwrite")
    return(gwas)
  }
  print("populating RSID...")
  gwas$chrbp <- paste0(gwas$CHR, ":", gwas$BP)

  chrpos_to_rsid <- vroom::vroom(marker_to_rsid_file, col_select=c(reference_build, "rsid"))
  gwas$RSID <- chrpos_to_rsid$rsid[match(gwas$chrbp, chrpos_to_rsid[[reference_build]])]
 
  return(gwas)
}


#neat faster way to read a big file: vroom(pipe("grep -w UA flights.tsv"), col_names = names(flights))
# or maybe need to pipe("awk 'NR == 1 || /Incoming/' foo.csv") to include the header row as well
read_and_format <-function(file_gwas, input_version=default, output_version="default", custom_file_columns=NA){
  if (input_version == "ieu_ukb_gwas_pipeline"){
    gwas <- vroom(file_gwas) %>%
      change_column_names(., ieu_ukb_gwas_columns) %>%
      standardise_columns(.)
  }
  else if (input_version == "metal") {
    gwas <- vroom(file_gwas) %>%
      change_column_names(., metal_columns) %>% 
      standardise_columns(.)
  }
  else if (input_version == "ukb_neale") {
    # UKB data produced by Neale lab
    print("reading variants")
    variants <- vroom(paste0("meta/neale_lab_variants_rsid_only.tsv")) # this is a col subset of variants.tsv.bgz from Neale lab
    print("reading gwas")
    gwas <-vroom(file_gwas, col_select = c("variant","minor_allele","minor_AF","beta","se","pval"))
    
    print("joining")
    if (nrow(gwas) == nrow(variants)) {
      merged <- left_join(gwas, variants, by = c("variant"="variant", "minor_allele"="alt")) %>%
        select(-variant) %>% select(SNP=rsid, everything())
    }
    
    print("formatting")
    gwas <- format_data(merged, type="outcome",
                        snp_col = "SNP",
                        beta_col = "beta",
                        se_col = "se",
                        effect_allele_col = "minor_allele",
                        other_allele_col = "ref",
                        eaf_col = "minor_AF",
                        pval_col = "pval")
    
    
  } else if (input_version == 'gwas_cat') {
    # harmonised data fromat from GWAS catalog
    gwas <- vroom(file_gwas,
                  col_select = c("variant_id",'beta', "standard_error","effect_allele","other_allele","effect_allele_frequency","p_value")) %>%
      format_data(., type="outcome",
                  snp_col = "variant_id",
                  beta_col = "beta",
                  se_col = "standard_error",
                  effect_allele_col = "effect_allele",
                  other_allele_col = "other_allele",
                  eaf_col = "effect_allele_frequency",
                  pval_col = "p_value")
  } else if (input_version == "custom") {
    custom_file_columns <- custom_file_columns %>% select_if(~all(!is.na(.)))
    stopifnot(dim(custom_file_columns)[2]>0)
    # data produced by IEU GWAS pipeline: pval col P_BOLT_LMM
    gwas <- vroom(file_gwas,
                  col_select = c(custom_file_columns[1,])) %>%
      format_data(.,
                  phenotype_col= custom_file_columns$trait,
                  snp_col = custom_file_columns$custom_SNP,
                  beta_col = custom_file_columns$custom_beta,
                  se_col = custom_file_columns$custom_se,
                  effect_allele_col = custom_file_columns$custom_effect_allele,
                  other_allele_col = custom_file_columns$custom_other_allele,
                  eaf_col = custom_file_columns$custom_eaf,
                  pval_col = custom_file_columns$custom_pval,
                  samplesize_col = custom_file_columns$custom_samplesize,
                  ncase_col = custom_file_columns$custom_ncase ,
                  ncontrol_col = custom_file_columns$custom_ncontrol,
                  chr_col = custom_file_columns$custom_chr,
                  pos_col = custom_file_columns$custom_pos,
                  z_col = custom_file_columns$custom_z,
                  info_col = custom_file_columns$custom_info,
                  units_col = custom_file_columns$custom_units,
                  gene_col = custom_file_columns$custom_gene
      )
    #TODO drop .outcome
  }
  
  return(gwas)
  
}
