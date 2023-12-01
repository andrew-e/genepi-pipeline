reference_builds = list(GRCh36="GRCh36", GRCh37="GRCh37", GRCh38="GRCh38")
available_liftover_conversions <- list()
available_liftover_conversions[[paste0(reference_builds$GRCh36, reference_builds$GRCh37)]]<-"hg18ToHg19.over.chain.gz"
available_liftover_conversions[[paste0(reference_builds$GRCh38, reference_builds$GRCh37)]]<-"hg38ToHg19.over.chain.gz"
available_liftover_conversions[[paste0(reference_builds$GRCh37, reference_builds$GRCh38)]]<-"hg19ToHg38.over.chain.gz"
available_liftover_conversions[[paste0(reference_builds$GRCh37, reference_builds$GRCh36)]]<-"hg19ToHg18.over.chain.gz"


convert_reference_build_via_liftover <- function(gwas, input_reference_build, output_reference_build, output_file) {
  gwas <- get_file_or_dataframe(gwas)
    
  if (!all(c(input_reference_build, output_reference_build) %in% reference_builds)) {
    stop("Error: invalid/unsupported reference builds when attempting to convert via liftover")
  }

  gwas$CHRBP <- paste(gwas$CHR, gwas$BP, sep=":")
  if (input_reference_build == reference_builds$GRCh37) {
    gwas$BP37 <- gwas$BP
  } else if (input_reference_build == reference_builds$GRCh38) {
    gwas$BP38 <- gwas$BP
  } else if (input_reference_build == reference_builds$GRCh36) {
    gwas$BP36 <- gwas$BP
  }

  bed_file_input <- tempfile(fileext = ".bed")
  bed_file_output <- tempfile(fileext = ".bed")
  unmapped <- tempfile(fileext = "unmapped")

  create_bed_file_from_gwas(gwas, bed_file_input)
  run_liftover(bed_file_input, bed_file_output, input_reference_build, output_reference_build, unmapped)
  gwas <- use_bed_file_to_update_gwas(gwas, bed_file_output)

  if(!missing(output_file) && shiny::isTruthy(output_file)) vroom::vroom_write(gwas, output_file)
  return(gwas)
}


create_bed_file_from_gwas <- function(gwas, output_file, marker_column="SNP") {
  gwas <- get_file_or_dataframe(gwas)

  bed_format <- tibble::tibble(
    CHR = paste0("chr", gwas$CHR),
    BP1 = gwas$BP,
    BP2 = gwas$BP+1,
    CHRBP = paste0(CHR, ":", BP1, "-", BP2)
  )

  vroom::vroom_write(bed_format, output_file, col_names=F, delim=" ")
  return(bed_format)
}


run_liftover <- function(bed_file_input, bed_file_output, input_build, output_build, unmapped) {
  lifover_binary <- paste0(liftover_dir, "liftOver")
  
  liftover_conversion <- available_liftover_conversions[[paste0(input_build, output_build)]]

  if (is.null(liftover_conversion)) {
    stop(paste("Error: liftOver format", output_format, "not recognized."))
  }

  chain_file <- paste0(liftover_dir, liftover_conversion)
  liftover_command <- paste(lifover_binary, bed_file_input, chain_file, bed_file_output, unmapped)
  print(liftover_command)
  system(liftover_command, wait=T)
}


use_bed_file_to_update_gwas <- function(gwas, bed_file) {
  liftover_bed <- vroom::vroom(bed_file, col_names=F)
  liftover_bed$X1 <- gsub("chr", "", liftover_bed$X1)
  liftover_bed$X4 <- sub("chr\\d:(\\d+)-\\d+", "\\1", liftover_bed$X4, perl=T)

  bed_map <- tibble::tibble(
    NEW_BP = liftover_bed$X2,
    ORIGINAL_CHRBP = paste(liftover_bed$X1, liftover_bed$X4, sep=":")
  )

  gwas <- merge(gwas, bed_map, by.x="CHRBP", by.y="ORIGINAL_CHRBP", all.x=T) |>
    dplyr::mutate(CHRBP=NULL, BP=NEW_BP)
  return(gwas)
}
