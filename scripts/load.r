#maybe look at changing this to creating a package: https://stackoverflow.com/a/48094346/6104011
source("../R/util.r")
source("../R/gwas_formatting.r")
source("../R/gwas_comparisons.r")
source("../R/collider_bias.r")
source("../R/graphs.r")
source("../R/coloc.r")
source("../R/mr.r")
source("../R/data_conversions.r")

genomic_data_dir <- if (Sys.getenv("GENOMIC_DATA_DIR") == "") "/mnt/storage/private/mrcieu/data/genomic_data/" else Sys.getenv("GENOMIC_DATA_DIR")
qtl_top_hits_dir  <- if (Sys.getenv("QTL_TOP_HITS_DIR") == "") "/mnt/storage/private/mrcieu/data/qtl_top_hits/" else Sys.getenv("QTL_TOP_HITS_DIR")

pqtl_top_hits_dir <- paste0(qtl_top_hits_dir, "/pqtl")
metabrain_top_hits_dir <- paste0(qtl_top_hits_dir, "/metabrain/top_hits")
metabrain_gwas_dir <- paste0(qtl_top_hits_dir, "/metabrain/gwas")
