source("load.r")

library(argparser, quietly = TRUE)

parser <- arg_parser("Create a markdown file of results")

parser <- add_argument(parser, "--markdown_file",
                       help = "Params to pass into markdown, a=b,c=d",
                       type = "character"
)
parser <- add_argument(parser, "--params",
                       help = "Params to pass into markdown, a=b,c=d",
                       type = "character"
)
output_file <- add_argument(parser, "--output_file",
                            help = "Name of output file",
                            type = "character"
)

args <- parse_args(parser)

rmd_params <- split_string_into_list(parser$params)

create_rmd_file(args$markdown_file, rmd_params, args$output_file)
