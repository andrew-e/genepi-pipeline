source("/home/scripts/functions/load.r")

library(argparser, quietly=TRUE)

p <- arg_parser("Create a manhattan plot from GWAS summary statistics")

p <- add_argument(p, "--gwas_filename", help="number to round", type="character")
p <- add_argument(p, "--manhattan_filename", help="number of decimal places", type="character")
p <- add_argument(p, "--qq_filename", help="include qq plot",  type="character")
p <- add_argument(p, "--include_qq", help="include qq plot", default=TRUE, type="logical")

args <- parse_args(p)

manhattan_and_qq(args$gwas_filename, args$manhattan_filename, args$qq_filename)
