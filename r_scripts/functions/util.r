DEFAULT_COLUMNS <- list(MARKER="MARKER", EA="EA", NONEA="NONEA", EAF="EAF", P="P", BETA="BETA", SE="SE")
METAL_COLUMNS <- list(MARKER="MarkerName", EA="Allele1", NONEA="Allele2", EAF="Freq1", P="P-value", BETA="Effect", SE="StdErr")
PLINK_COLUMNS = list(MARKER="ID", EA="REF", NONEA="ALT1", EAF="Freq1", P="P.value", EFFECT="BETA", SE="StdErr", CHR="CHROME", BP="POS")

dbsnp.hg37 <- "/mnt/storage/private/mrcieu/research/mr-eve/vcf-reference-datasets/dbsnp/dbsnp.v153.b37.db"

save_file <- function(df, filename, save_location, format) {

}

chrpos_to_rsid <- function(gwas_filename, column = "SNP") {
}

rsid_to_chrpos <- function(gwas_filename, column = "SNP") {
    gwas <- data.table::fread(gwas_filename)
    gwas[[column]] <- gsub('^rs', '', gwas[[column]])

    rsid_list <- toString(gwas[[column]])

    sqlite_query <- paste0("SELECT 'rs' || rsid || ',' || chrom || ':' || coord FROM rsid_to_coord WHERE rsid IN (", rsid_list, ")")

    sqlite_command <- paste("sqlite3", dbsnp.hg37, sqlite_query)
    output <- system(sqlite_command, intern = TRUE)
}


column_names <- function(columns = list()) {
    chosen_names <- DEFAULT_COLUMNS

    for (column in columns) {
        name <- columns[[column]]
        chosen_names[[name]] <- if(!is.null(columns[[name]])) columns[[name]]
    }
    return(chosen_names)
}

#' change_header_names: change header names of file
#'
#' @param gwas_filename: a dataframe that includes CHR, CP, P, and SNP
#' @param column_names: named list of columns to change (ex. MARKER="ID")
#' @return: changes the gwas_filename inline.
#' @examples
change_header_names <- function(gwas_filename, column_names = list()) {
    if (length(column_names) != 0) {
        for (name in names(column_names)) {
            sed_command <- paste0("sed -i '0,/", column_names[[name]],
                                  "/{s/", column_names[[name]], "/", name, "/}' ",
                                  gwas_filename)
            system(sed_command)
        }
    }
}


#' standardise_gwas_data: standardise gwas column names and columns
#'
#' @param gwas_filename: gwas file to read in
#' @param sep: separator for gwas file
#' @param column_names: names of columns that need to be renamed
#'    (can use helper lists at the top of the file (METAL_COLUMNS)
#' @param separate_chr_bp: will split MARKER into separate CHR and BP columns
#' @param save_as: file name and location to save to (will update to
#' @return 2 plots: one manhattan plot and one QQ plot (with lambda included)
#' @examples
standardise_gwas_data <- function(gwas_filename=gwas_filename, sep="\t", column_names=list(), separate_chr_bp=F, save_as=F) {
    change_header_names(gwas_filename = gwas_filename, column_names = column_names)

    #example of changing  chr:pos:ea:oa to leave only chr:pos
    #r <- regexpr("([0-9]+:[0-9]+)*", pqtl_data$SNPID)
    #pqtl_data$SNPID_STRIPPED <- NA
    #pqtl_data$SNPID_STRIPPED[which(r != -1)] <- regmatches(pqtl_data$SNPID, r)

    gwas_data <- read.table(gwas_filename, sep=sep, header=T)

    if (separate_chr_bp) {
        gwas_data <- separate(data = gwas_data, col = "MARKER", into = c("CHR", "BP"), sep = ":", remove = F)
    }

    if (join_chr_bp) {
        gwas_data <- separate(data = gwas_data, col = "MARKER", into = c("CHR", "BP"), sep = ":", remove = F)
    }

    if (isTruthy(save_as)) {
        write.table(gwas_data, save_as, row.names = F, quote = F)
    } else {
        write.table(gwas_data, gwas_filename, row.names = F, quote = F)
    }
}
