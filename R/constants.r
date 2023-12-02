number_of_cpus_available <- as.numeric(get_env_var("SLURM_CPUS_ON_NODE", 1))
genomic_data_dir <- get_env_var("GENOMIC_DATA_DIR", "/mnt/storage/private/mrcieu/data/genomic_data/")
qtl_top_hits_dir <- get_env_var("QTL_TOP_HITS_DIR", "/mnt/storage/private/mrcieu/data/qtl_top_hits/")
liftover_dir <- get_env_var("LIFTOVER_DIR", "/projects/MRC-IEU/programs/twosampleMR/liftover/")

pqtl_top_hits_dir <- paste0(qtl_top_hits_dir, "/pqtl")
metabrain_top_hits_dir <- paste0(qtl_top_hits_dir, "/metabrain/top_hits")
metabrain_gwas_dir <- paste0(qtl_top_hits_dir, "/metabrain/gwas")
