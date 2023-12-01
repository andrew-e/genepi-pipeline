#!/bin/bash

#SBATCH --account=<smed>
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

singularity run -B $(pwd)/R:/home/R -B $(pwd)/scripts:/home/scripts -B /user/work/$(whoami) -B /user/home/$(whoami) -B /mnt/storage/private/mrcieu/data/ --env-file .env --pwd /home/scripts docker://andrewrrelmore/genepi_pipeline:refactor $SCRIPT
