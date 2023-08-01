dbsnp.hg37 <- "/mnt/storage/private/mrcieu/research/mr-eve/vcf-reference-datasets/dbsnp/dbsnp.v153.b37.db"

split_string_into_vector <- function(input_string) {
  return(unlist(strsplit(input_string, ',')))
}

#' vroom_region: faster way of opening a GWAS, only load a specific chromosome's worth of data
#'
#' NOTE: only works with data that has been standardised, through `standardise_gwas`, or at least a tsv
vroom_chr <- function(gwas_file, chr, col_select=NULL) {
  pipe_command <- paste0("head -n1 ", gwas_file, " && grep '\t", chr, "\t' ", gwas_file)

  gwas <- vroom::vroom(pipe(pipe_command, col_select = col_select))
  return(gwas)
}

create_dir_for_files <- function(...) {
  filenames <- list(...)
  library(stringr)

  for (filename in filenames) {
    filepath <- file.path(str_extract(filename, "^(.*/)"))
    if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)
  }
}

map_rsid_list_to_snps <- function(gwas, rsids=c()) {
  gwas <- subset(gwas, RSID %in% rsids)
  return(gwas$SNP)
}

#TODO: unused, buy maybe helpful in the future
rsid_to_chrpos <- function(gwas_filename, column = "SNP") {
  gwas <- data.table::fread(gwas_filename)
  gwas[[column]] <- gsub('^rs', '', gwas[[column]])

  rsid_list <- toString(gwas[[column]])

  sqlite_query <- paste0("SELECT 'rs' || rsid || ',' || chrom || ':' || coord FROM rsid_to_coord WHERE rsid IN (", rsid_list, ")")

  sqlite_command <- paste("sqlite3", dbsnp.hg37, sqlite_query)
  output <- system(sqlite_command, intern = TRUE)
}

