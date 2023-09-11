#maybe look at changing this to creating a package: https://stackoverflow.com/a/48094346/6104011

source("functions/util.r")
source("functions/gwas_formatting.r")
source("functions/gwas_comparisons.r")
source("functions/collider_bias.r")
source("functions/graphs.r")
source("functions/coloc.r")

#TODO: change these dir references to populate from environment variables, and fallback to some default
thousand_genomes_dir <- "/user/work/wt23152/genome_data/1000genomes/"
pqtl_top_hits_dir <- "/user/work/wt23152/pqtl_data/"