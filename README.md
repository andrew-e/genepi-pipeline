# genepi-pipeline

Work in progress: a pipeline for creating reusable steps for genetic epideiology projects at the Univsersity of Bristol.

Goals:
* Reproducible reserach
* Reusable pipelines that can be used on different projects
* Shared code and steps that can be updated according to latest knowledge and practices

## Onboarding

Steps to start using the pipeline
#### Clone the repository on blue cristal 4
`git clone git@github.com:andrew-e/genepi-pipeline.git && cd genepi-pipeline`

#### Create and activate conda environment
`conda env create --file environment.yml`
`conda activate genepi-pipeline`

#### Populate the your DATA_DIR and RESULTS_DIR environment variables in .env file
These should probably be in your *work* or *scratch* space (`/user/work$(whoami)/...`)

These can't be your RDFS space, as that volume won't mount to into the singularity container

#### Run the test pipeline
`snakemake --profile snakemake/`

## How it works:

There are 2 main components to the pipeline
1. Snakemake to define the steps to complete for each pipeline.  Each snakemake step spins up
2. Docker / Singularity container with intalled (and version pinned) languages (R and python), packages, os libraries, and code
3. Slurm: each snakemake step spins up a singularity container inside a slurm job.  Each step can specify different resource requirements

We can define different Snakefiles for different pipelines as it becomes more mature


## Repository Organisation

* The `Dockerfile` and `docker/` directory hold the information for creating the docker image that the pipeline runs
* `r_scripts` code that can be called from snakemake
    * `r_scripts/functions/` holds r code that can be called and reused by any step in the pipeline (accessed by a cli script)
    * `r_scripts/cli/` holds the scripts that can be easily called by snakemake (`Rscript example.r --input_file example.txt`)
