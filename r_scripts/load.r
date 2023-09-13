#maybe look at changing this to creating a package: https://stackoverflow.com/a/48094346/6104011
source("functions/util.r")
source("functions/gwas_formatting.r")
source("functions/gwas_comparisons.r")
source("functions/collider_bias.r")
source("functions/graphs.r")
source("functions/coloc.r")

genome_data_dir <- if (is.na(Sys.getenv("THOUSAND_GENOMES_DIR"))) "/mnt/storage/private/mrcieu/data/1000genomes/" else Sys.getenv("THOUSAND_GENOMES_DIR")
qtl_top_hits_dir  <- if (is.na(Sys.getenv("QTL_TOP_HITS_DIR"))) "/mnt/storage/private/mrcieu/data/qtl_top_hits/" else Sys.getenv("QTL_TOP_HITS_DIR")