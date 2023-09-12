#!/bin/bash
set -e

conda activate genepi-pipeline
module load apps/singularity/3.8.3
tmux new -s genepi-pipeline
