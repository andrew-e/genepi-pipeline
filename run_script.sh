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

singularity run -B /user/work/$USER -B /user/home/$USER -pwd /home/r_scripts/ docker://andrewrrelmore/genepi_pipeline:latest $SCRIPT