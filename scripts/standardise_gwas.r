source("load.r")
library(argparser, quietly = TRUE)

parser <- arg_parser("Standardise a GWAS to be ready for the rest of the pipeline")

parser <- add_argument(parser, "--input_gwas",
                       help = "Comma separated list of filenames of GWASes to standardise",
                       type = "character"
)
parser <- add_argument(parser, "--input_columns",
                       help = "Map of column names for pipeline to change it to",
                       type = "character"
)
parser <- add_argument(parser, "--input_format",
                       help = "Input format of the gwas (ie. metal, bolt, plink, default)",
                       type = "character",
                       default = "default"
)
parser <- add_argument(parser, "--output_gwas",
                       help = "Comma separated list of filenames of the standardised GWASes",
                       type = "character"
)
parser <- add_argument(parser, "--populate_rsid",
                       help = "Should GWAS Populate RSID (based on 1000 genomes data)",
                       type = "character",
                       default = "PARTIAL"
)
parser <- add_argument(parser, "--to-output",
                       help = "Flag to format the standardised GWAS into a different output",
                       flag = T
)

args <- parse_args(parser)
create_dir_for_files(args$output_gwas)

if (!args$to_output) {
  bespoke_column_map <-split_string_into_named_list(args$input_columns)
  standardise_gwas(args$input_gwas,
                   args$output_gwas,
                   input_format = args$input_format,
                   populate_rsid_option = args$populate_rsid,
                   bespoke_column_map = bespoke_column_map
  )
} else {
    format_gwas_output(args$input_gwas, args$output_gwas, args$input_format)
    gc()
}
