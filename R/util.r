get_file_or_dataframe <- function(input, columns=NULL, snps=NULL) {
  if (is.data.frame(input)) {
    input <- dplyr::select(input, `if`(is.null(columns), dplyr::all_of(dplyr::everything()), dplyr::all_of(columns))) |>
        subset(`if`(is.null(snps), T, SNP %in% snps))
  }
  else {
    if (!file.exists(input)) stop(paste("Error:", input, "can't be found"))
    if (!is.null(snps)) {
      input <- vroom_snps(input, snps) |>
          dplyr::select(`if`(is.null(columns), dplyr::all_of(dplyr::everything()), dplyr::all_of(columns)))
    }
    else {
      if (is.null(columns)) {
        input <- vroom::vroom(input)
      } else {
        input <- vroom::vroom(input, col_select = dplyr::all_of(columns))
      }
    }
    input <- subset(input, `if`(is.null(snps), T, SNP %in% snps))
  }
  return(input)
}

#' vroom_snps: If you only need to get a handful of SNPs out of a whole GWAS, this saves time and memory
#' NOTE: only works with data that has been standardised, through `standardise_gwas`, or at least a tsv
vroom_snps <- function(gwas_file, snps=c()){
  snps <- paste(snps, collapse="\t|")

  if (endsWith(gwas_file, ".gz")) {
    if (Sys.info()["sysname"] == "Darwin") {
      grep_command <- paste0("zcat < ", gwas_file, " | head -n 1 && rg -Iz '", snps, "' ", gwas_file)
    } else {
      grep_command <- paste0("zcat ", gwas_file, " | head -n 1 && rg -Iz '", snps, "' ", gwas_file)
    }
  }
  else {
    grep_command <- paste0("head -n 1", gwas_file, " && rg -I '", snps, "' ", gwas_file)
  }

  snps_in_gwas <- vroom::vroom(pipe(grep_command))
  return(snps_in_gwas)
}

gwas_region <- function(gwas, chr, bp, range = 500000) {
  return(dplyr::filter(gwas, CHR == chr &BP > (bp - floor(range/2)) & BP < (bp + floor(range/2))))
}


split_string_into_vector <- function(input_string) {
  return(unlist(strsplit(input_string, '[ ]')))
}

parse_gwas_input_column_maps <- function(input_column_string) {
  column_map_as_a_string <- unlist(strsplit(input_column_string, '[:]'))
}

extract_numbers_from_string <- function(string) {
  regmatches(string, gregexpr("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?", string))
}

split_string_into_named_list <- function(input_string) {
  split <- unlist(strsplit(input_string, '[=,]'))
  names <- split[c(T, F)]
  values <- split[c(F, T)]

  return(structure(as.list(values), names=names))
}

get_env_var <- function(env_var_name, default_value) {
  if (Sys.getenv(env_var_name) == "") {
    return(default_value)
  }
  else {
    return(Sys.getenv(env_var_name))
  }
}

file_prefix <- function(file_path) {
  file_name <- basename(file_path)
  file_prefix <- sub("\\..*", "", file_name)
  file_prefix <- sub("_std", "", file_prefix)
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
  temp_file <- tempfile(fileext = ".rmd", tmpdir = "/tmp")
  file.copy(rmd_file, temp_file, overwrite = TRUE)

  rmarkdown::render(temp_file,
                    output_file = output_file,
                    params = params,
                    quiet = T
  )
}

get_other_docker_tag <- function() {
  tag_to_match <- "test"
  docker_url <- "https://hub.docker.com/v2/repositories/andrewrrelmore/genepi_pipeline/tags/"
  tag_information <- c()
  tryCatch(
   expr = {
      response <- httr::GET(docker_url, httr::accept_json())
      tag_information <- httr::content(response, type="application/json")$results
   },
   error = function(e) {
     message(paste("Call to docker hub failed:", e))
     return(tag_to_match)
   }
  )

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
