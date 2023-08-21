source("load.r")
library(argparser, quietly = TRUE)

parser <- arg_parser("Comparison of expected vs observed replication rates")

parser <- add_argument(parser,
                       "--gwas_filenames",
                       help = "Comma separated list of GWAS filenames",
                       type = "character"
)
parser <- add_argument(parser,
                       "--clumped_filenames",
                       help = "Comma separated list of clumped SNP filenames",
                       type = "character"
)
parser <- add_argument(parser,
                       "--result_output",
                       help = "Output file name of results from gwas comparison",
                       type = "character"
)
parser <- add_argument(parser,
                       "--variants_output",
                       help = "Output file name of variants from gas comparison",
                       type = "character"
)

args <- parse_args(parser)

gwas_filenames <- split_string_into_vector(args$gwas_filenames)
clumped_filenames <- split_string_into_vector(args$clumped_filenames)

if (length(gwas_filenames) != length(clumped_filenames)) {
  stop("List size of --gwas_filenames and --clumped_filenames are different.  They must be the same.")
}

compare_replication_across_all_gwas_permutations(gwas_filename, clumped_filenames, args$result_output, args$variants_output)
