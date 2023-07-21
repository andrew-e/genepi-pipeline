source("home/r_scripts/load.r")
library(argparser, quietly = TRUE)

parser <- arg_parser("Correct for Collider Bias between incidence and subsequent GWASes")

parser <- add_argument(parser, "--incidence_gwas",
                       help = "Indience GWAS",
                       type = "character"

)
parser <- add_argument(parser, "--subsequent_gwas",
                       help = "Subsequent GWAS",
                       type = "character"
)
parser <- add_argument(parser, "--clumped_file",
                       help = "Clumped SNP list to run collider bias corrections against",
                       type = "character"
)
parser <- add_argument(parser, "--collider_bias_results_output",
                       help = "Output file of Collider Bias Slopes", 
                       type = "character"
)
parser <- add_argument(parser, "--collider_bias_adjusted_output",
                       help = "Output file of Collider Bias corrected GWAS results",
                       type = "character"
)
parser <- add_argument(parser, "--p_value_threshold", 
                       help = "p value threshold to run corrections",
                       type = "numeric", default = 0.001
)
parser <- add_argument(parser, "--include_slopehunter",
                       help = "Run slopehunter collider bias correction",
                       type = "logical", default = T
)
parser <- add_argument(parser, "--include_dudbridge",
                       help = "Run dudbridge collider bias correction",
                       type = "logical", default = T
)

args <- parse_args(parser)
p_value_threshold <- as.numeric(args.p_value_threshold)

correct_for_collider_bias(args.incidence_gwas,
                          args.subsequent_gwas,
                          args.clumped_file,
                          args.collider_bias_results_output,
                          args.collider_bias_adjusted_output,
                          p_value_threshold,
                          args.include_slopehunter,
                          args.include_dudbridge
)
