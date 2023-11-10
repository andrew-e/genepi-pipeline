#!/bin/bash
set -e
SMK_FILE=$1

conda activate genepi-pipeline
module load apps/singularity/3.8.3
tmux new-session -d -s my_session "snakemake --snakefile ${SMK_FILE} --profile snakemake/bc4/"