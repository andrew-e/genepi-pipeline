dbsnp.hg37 <- "/mnt/storage/private/mrcieu/research/mr-eve/vcf-reference-datasets/dbsnp/dbsnp.v153.b37.db"

split_string_into_vector <- function(input_string) {
  return(unlist(strsplit(input_string, ', ')))
}

split_string_into_list <- function(input_string) {
  return(unlist(strsplit(input_string, '=,')))
}

#' vroom_chr: faster way of opening a GWAS, only load a specific chromosome's worth of data
#'
#' NOTE: only works with data that has been standardised, through `standardise_gwas`, or at least a tsv
vroom_chr <- function(gwas_file, chr, col_select=NULL) {
  pipe_command <- paste0("head -n1 ", gwas_file, " && rg -Iz '\t", chr, "\t' ", gwas_file)

  gwas <- vroom::vroom(pipe(pipe_command, col_select = col_select))
  return(gwas)
}

gwas_region <- function(gwas, chr, bp, range = 250000) {
  return(subset(gwas, CHR == chr & BP > (bp - range) & BP < (bp + range)))
}

file_prefix <- function(file_path) {
  file_name <- basename(file_path)
  file_prefix <- sub("\\..*", "", file_name)
  return(file_prefix)
}

create_dir_for_files <- function(...) {
  filenames <- list(...)
  library(stringr, quietly = T)

  for (filename in filenames) {
    filepath <- file.path(stringr::str_extract(filename, "^(.*/)"))
    if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)
  }
}

map_rsid_list_to_snps <- function(gwas, rsids=c()) {
  gwas <- subset(gwas, RSID %in% rsids)
  return(gwas$SNP)
}

create_rmd_file <- function(rmd_file, params = list(), output_file) {
  library(rmarkdown, quietly = T)
  rmarkdown::render(paste0("markdown/", rmd_file),
                    "pdf_document",
                    output_file = output_file,
                    params = params)
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
