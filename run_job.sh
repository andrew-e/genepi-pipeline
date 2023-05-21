#!/bin/bash

#SBATCH --account=smed001801
#SBATCH --partition=mrcieu
#SBATCH --job-name=singularity_job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --output=/user/work/%u/slurm/pipeline.out

module add apps/singularity/3.8.3

SCRATCH_MNT_DIR=/user/work/$USER
SCRIPT=""

if [ "$SCRIPT" ]
then
    SCRIPT="Rscript $SCRIPT"
fi

singularity run -B $SCRATCH_MNT_DIR:./scratch docker://andrewrrelmore/genepi_pipeline:latest $SCRIPT

