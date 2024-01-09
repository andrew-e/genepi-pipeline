#!/bin/bash

#SBATCH --account=smed001801
#SBATCH --partition=mrcieu
#SBATCH --job-name=bespoke_pipeline_job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --output=/user/work/%u/slurm_logs/bespoke_pipeline_job_%j.out

module add apps/singularity/3.8.3
USER=$(whoami)
#FILES="/user/work/wt23152/005/data/gwas/scz_gwas_filtered_all_std.tsv.gz,/user/work/wt23152/005/data/gwas/ebi-a-GCST003156_std.tsv.gz,/user/work/wt23152/005/data/gwas/ieu-b-30_std.tsv.gz,/user/work/wt23152/005/data/gwas/ieu-b-32_std.tsv.gz,/user/work/wt23152/005/data/gwas/ieu-b-33_std.tsv.gz,/user/work/wt23152/005/data/gwas/ieu-b-34_std.tsv.gz,/user/work/wt23152/005/data/gwas/ukb-a-248_std.tsv.gz,/user/work/wt23152/005/data/gwas/ukb-a-448_std.tsv.gz,/user/work/wt23152/005/data/gwas/ukb-b-12039_std.tsv.gz,/user/work/wt23152/005/data/gwas/ukb-b-4226_std.tsv.gz"
NS="150000,14267,36410,35266,36242,35212,336107,180203,454893,463010"
FILES="/user/work/wt23152/005/data/gwas/scz_gwas_filtered_all_std.tsv.gz,/user/work/wt23152/005/data/gwas/ieu-b-30_std.tsv.gz,/user/work/wt23152/005/data/gwas/ieu-b-32_std.tsv.gz,/user/work/wt23152/005/data/gwas/ieu-b-33_std.tsv.gz,/user/work/wt23152/005/data/gwas/ieu-b-34_std.tsv.gz,/user/work/wt23152/005/data/gwas/ukb-a-248_std.tsv.gz,/user/work/wt23152/005/data/gwas/ukb-a-448_std.tsv.gz,/user/work/wt23152/005/data/gwas/ukb-b-12039_std.tsv.gz,/user/work/wt23152/005/data/gwas/ukb-b-4226_std.tsv.gz"
NS="150000,36410,35266,36242,35212,336107,180203,454893,463010"
SCRIPT="./run_ldsc.sh $FILES $NS EUR /user/work/wt23152/005/results/ldsc/result_eur"

singularity run -B $(pwd)/R:/home/R -B $(pwd)/scripts:/home/scripts -B /user/work/$(whoami) -B /user/home/$(whoami) -B /mnt/storage/private/mrcieu/data/ --env-file .env --pwd /home/scripts docker://andrewrrelmore/genepi_pipeline:test $SCRIPT
