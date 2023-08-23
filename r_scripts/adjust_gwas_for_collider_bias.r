source("load.r")
library(argparser, quietly = TRUE)

parser <- arg_parser("Correct for Collider Bias between incidence and subsequent GWASes")

parser <- add_argument(parser, "--gwas",
                       help = "GWAS to be adjsuted by collider bias beta and se",
                       type = "character"
)
parser <- add_argument(parser, "--harmonised_effects_gwas",
                       help = "Harmonised Effects GWAS, created by SlopeHunter",
                       type = "character"
)
parser <- add_argument(parser, "--beta",
                       help = "beta value to adjusted the GWAS by",
                       type = "character"
)
parser <- add_argument(parser, "--se",
                       help = "Clumped SNP list to run collider bias corrections against",
                       type = "character"
)
parser <- add_argument(parser, "--output_gwas",
                       help = "Output GWAS file of collider bias adjusted result",
                       type = "character"
)

args <- parse_args(parser)

type <- "bespoke_adjustment"
adjust_gwas_data_from_weights_and_save(args$harmonised_effects_gwas,
                                       type,
                                       args$beta,
                                       args$se,
                                       args$output_gwas
)
