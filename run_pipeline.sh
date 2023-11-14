#!/bin/bash
set -e
SMK_FILE=$1

#if [$(hostname)...]

conda activate genepi-pipeline || echo "Error: Please ensure your conda environment is initialised and 'genepi-pipeline' env is available" && exit 1
module load apps/singularity/3.8.3
tmux new-session -d -s my_session "snakemake --snakefile ${SMK_FILE} --profile snakemake/bc4/"