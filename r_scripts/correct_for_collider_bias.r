source("home/r_scripts/load.r")
library(argparser, quietly = TRUE)

parser <- arg_parser("Correct for Collider Bias between incidence and subsequent GWASes")
parser <- add_argument(parser, "--incidence_gwas", help = "Indience GWAS", type = "character")
parser <- add_argument(parser, "--subsequent_gwas", help = "Subsequent GWAS", type = "character")
parser <- add_argument(parser, "--clumped_snps", help = "Clumped SNP list to run collider bias corrections against", type = "character")
parser <- add_argument(parser, "--output_file", help = "Output file of Collider Bias corrected GWAS results", type = "character")
parser <- add_argument(parser, "--p_value_thresholds", help = "Comma separated list of p value thresholds to run corrections", type = "character", default = 0.001)
parser <- add_argument(parser, "--include_slopehunter", help = "Run slopehunter collider bias correction", type = "logical", default = T)
parser <- add_argument(parser, "--include_dudbridge", help = "Run dudbridge collider bias correction", type = "logical", default = T)

args <- parse_args(parser)

p_value_thresholds <- as.numeric(split_string_into_list(args.p_value_thresholds))
correct_for_collider_bias(args.incidence_gwas, args.subsequent_gwas, args.clumped_snps, args.output_file, p_value_thresholds, args.include_slopehunter, args.include_dudbridge)
