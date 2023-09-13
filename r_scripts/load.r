#maybe look at changing this to creating a package: https://stackoverflow.com/a/48094346/6104011
source("functions/util.r")
source("functions/gwas_formatting.r")
source("functions/gwas_comparisons.r")
source("functions/collider_bias.r")
source("functions/graphs.r")
source("functions/coloc.r")

genome_data_dir <- if (is.na(Sys.getenv("1000_GENOMES_DIR"))) "/mnt/storage/private/mrcieu/data/1000genomes/" else Sys.getenv("1000_GENOMES_DIR")
pqtl_top_hits_dir <- "/user/work/wt23152/pqtl_data/"