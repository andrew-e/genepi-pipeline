#!/bin/bash
set -e

#if conda info --envs | grep -q "genepi-pipeline"; then echo ""; else conda create environment.yaml; fi
if conda info --envs | grep -q "genepi-pipeline"; then echo "conda environment ready"; else echo "install shittt"; fi
conda activate genepi-pipeline
module load apps/singularity/3.8.3
tmux new -s genepi-pipeline
