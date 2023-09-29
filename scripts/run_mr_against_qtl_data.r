source("load.r")
library(argparser, quietly = TRUE)

parser <- arg_parser("Comparison of expected vs observed replication rates")

parser <- add_argument(parser,
                       "--gwas_filenames",
                       help = "Comma separated list of GWAS filenames",
                       type = "character",
                       nargs = Inf
)
parser <- add_argument(parser,
                       "--ancestries",
                       help = "Comma separated list of ancestries related to GWAS",
                       type = "character",
                       nargs = Inf
)
parser <- add_argument(parser,
                       "--dataset",
                       help = "Dataset of type of QTL analysis to run ('pqtl', and 'metabrain' accepted)",
                       type = "character"
)
parser <- add_argument(parser,
                       "--result_output",
                       help = "Output file name of results from gwas comparison",
                       type = "character"
)

args <- parse_args(parser)
mr_function <- list(metabrain = perform_mr_on_metabrain_datasets, pqtl = perform_mr_on_pqtl_datasets)

gwas_filenames <- split_string_into_vector(args$gwas_filenames)
ancestries <- split_string_into_vector(args$ancestries)
result_outputs <- split_string_into_vector(args$result_outputs)

if (length(gwas_filenames) != length(ancestries) || length(gwas_filenames) != length(result_outputs)) {
  stop("List size of --gwas_filenames --ancestries and --result_outputs are different.  They must be the same.")
}

for (i in seq_along(input_gwases)) {
  create_dir_for_files(output_gwases[i])
  mr_function[[args$dataset]](input_gwases[i], ancestries[i], result_outputs[i])
}