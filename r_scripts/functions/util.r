split_string_into_vector <- function(input_string) {
  return(unlist(strsplit(input_string, '[ ]')))
}

parse_gwas_input_column_maps <- function(input_column_string) {
  column_map_as_a_string <- unlist(strsplit(input_column_string, '[:]'))
}

split_string_into_named_list <- function(input_string) {
  split <- unlist(strsplit(input_string, '[=,]'))
  names <- split[c(T, F)]
  values <- split[c(F, T)]

  return(structure(as.list(values), names=names))
}

#' vroom_chr: faster way of opening a GWAS, only load a specific chromosome's worth of data
#'
#' NOTE: only works with data that has been standardised, through `standardise_gwas`, or at least a tsv
vroom_chr <- function(gwas_file, chr, col_select=NULL) {
  pipe_command <- paste0("head -n1 ", gwas_file, " && rg -Iz '\t", chr, "\t' ", gwas_file)
  pipe_command <- paste0("rg -Iz '\t", chr, "\t' ", gwas_file)

  gwas <- vroom::vroom(pipe(pipe_command), col_select = col_select)
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

create_html_from_rmd <- function(rmd_file, params = list(), output_file) {
  library(rmarkdown, quietly = T)
  temp_file <- tempfile(fileext = ".Rmd")
  file.copy(rmd_file, temp_file, overwrite = TRUE)

  rmarkdown::render(temp_file,
                    output_file = output_file,
                    params = params)
}

get_other_docker_tag <- function() {
  tag_to_match <- "test"
  docker_url <- "https://hub.docker.com/v2/repositories/andrewrrelmore/genepi_pipeline/tags/"

  response <- httr::GET(docker_url, httr::accept_json())
  tag_information <- httr::content(response, type="application/json")$results

  #we can't be sure about tag order, so iterating over it twice
  for (tag in tag_information) {
    if (tag$name == tag_to_match) {
      digest <- tag$digest
    }
  }
  for (tag in tag_information) {
    if (tag$digest == digest & tag$name != tag_to_match) {
      return(tag$name)
    }
  }
  return(tag_to_match)
}

#TODO: unused, buy maybe helpful in the future
rsid_to_chrpos <- function(gwas_filename, column = "SNP") {
  dbsnp.hg37 <- "/mnt/storage/private/mrcieu/research/mr-eve/vcf-reference-datasets/dbsnp/dbsnp.v153.b37.db"
  gwas <- data.table::fread(gwas_filename)
  gwas[[column]] <- gsub('^rs', '', gwas[[column]])

  rsid_list <- toString(gwas[[column]])

  sqlite_query <- paste0("SELECT 'rs' || rsid || ',' || chrom || ':' || coord FROM rsid_to_coord WHERE rsid IN (", rsid_list, ")")

  sqlite_command <- paste("sqlite3", dbsnp.hg37, sqlite_query)
  output <- system(sqlite_command, intern = TRUE)
}
