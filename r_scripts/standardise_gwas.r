source("load.r")
library(argparser, quietly = TRUE)

parser <- arg_parser("Standardise a GWAS to be ready for the rest of the pipeline")

parser <- add_argument(parser, "--input_gwas",
                       help = "Comma separated list of filenames of GWASes to standardise",
                       type = "character",
                       nargs = Inf
)
parser <- add_argument(parser, "--input_columns",
                       help = "Map of column names for pipeline to change it to",
                       type = "character",
                       nargs = Inf
)
parser <- add_argument(parser, "--input_format",
                       help = "Input format of the gwas (ie. METAL, BOLT, plink)",
                       type = "character",
                       default = "default"
)
parser <- add_argument(parser, "--output_gwas",
                       help = "Comma separated list of filenames of the standardised GWASes",
                       type = "character",
                       nargs = Inf
)
parser <- add_argument(parser, "--populate_rsid",
                       help = "Should GWAS Populate RSID (based on 1000 genomes data)",
                       flag = T
)
parser <- add_argument(parser, "--to-output",
                       help = "Flag to format the standardised GWAS into a different output",
                       flag = T
)

args <- parse_args(parser)
input_gwases <- split_string_into_vector(args$input_gwas)
input_columns <- split_string_into_vector(args$input_columns)
output_gwases <- split_string_into_vector(args$output_gwas)

if (length(input_gwases) != length(output_gwases)) {
  stop("input_gwas and output_gwas need to be the same length")
}

if (!args$to_output) {
  for (i in seq_along(input_gwases)) {
    create_dir_for_files(output_gwases[i])
    standardise_gwas(input_gwases[i], output_gwases[i], args$input_format, args$populate_rsid, bespoke_column_map = input_columns[i])
    gc()
  }
} else {
  for (i in seq_along(input_gwases)) {
    create_dir_for_files(output_gwases[i])
    format_gwas_output(input_gwases[i], output_gwases[i], args$input_format)
    gc()
  }
}
