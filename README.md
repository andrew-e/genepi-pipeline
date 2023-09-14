# genepi-pipeline

Work in progress: a pipeline for creating reusable steps for genetic epidemiology projects at the University of Bristol.

Goals:
* Reproducible research
* Reusable pipelines that can be used on different projects
* Shared code and steps that can be updated according to the latest knowledge and practices

## Onboarding

Steps to start using the pipeline
### 1. Clone the repository on BlueCrystal4
`git clone git@github.com:andrew-e/genepi-pipeline.git && cd genepi-pipeline`

### 2. Create and activate conda environment
Ensure you [have conda installed](https://www.acrc.bris.ac.uk/protected/hpc-docs/software/python_conda.html) and initialised before running these commands:
```
conda env create --file environment.yml
conda activate genepi-pipeline
module load apps/singularity/3.8.3`
```
### 3. Populate .env and input.json files

### `cp .env_example .env`
* populate the DATA_DIR, RESULTS_DIR and RDFS_DIR environment variables in .env file
These should probably be in your *work* or *scratch* space (`/user/work/$(whoami)/...`)
* RDFS_DIR is optional, but very useful.  All generated files can be copied automatically.  Please ensure the path
ends in `working/`

### `cp workflow/input_<specific_pipeline>.json input.json`
* Each pipeline (as defined in `workflow` directory) has its own input format.  Please copy and populate the json file of the pipline you want to run
* With each GWAS, you can specify header names ex. `{"P":"your_gwas_pval_col", ...}`, if you do not specify header names it will assume your GWAS has the headers below.

### 4. Run the pipeline
### `snakemake --snakefile workflow/<specific_pipeline>.cmk --profile snakemake/`

If there are errors while running the pipeline, you can find error messages either directly on the screen, or there may be a slurm log file that is outputted on error, which you can look at for clues.

## Usage:

The standard column naming for GWASes are:

If you do not specify the expected format of the GWAS, all pipelines will assume the following headers:

|           | SNP | CHR | BP  | BETA | SE  | P   | EA  | OA  | EAF | RSID | OR  | OR_LB | OR_UB |
|-----------|-----|-----|-----|------|-----|-----|-----|-----|-----|:-----|-----|-------|-------|
| Mandatory |     | x   | x   | x    | x   | x   | x   | x   |     |      |     |       |       |

### GWAS Standardisation

All pipelines will standardise each GWAS before running the subsequent steps.  The `SNP` field will be recalculated as

`CHR:POS_EA_OA`, where EA and OA are ordered alphabetically, and the subsequent BETA and EAF will be adjusted accordingly


## How it works:

There are 3 main components to the pipeline
1. Snakemake to define the steps to complete for each pipeline.
2. Docker / Singularity container with installed languages (R and python), packages, os libraries, and code
3. Slurm: each snakemake step spins up a singularity container inside a slurm job.  Each step can specify different slurm requirements.

## Repository Organisation

* The `Dockerfile` and `docker/` directory hold the information for creating the docker image that the pipeline runs
* `r_scripts` code that can be called from the pipelines, inside the container
  * `r_scripts/` holds the scripts that can be easily called by snakemake (`Rscript example.r --input_file example.txt`)
  * `r_scripts/functions/` holds r code that can be called and reused by any step in the pipeline (accessed by a cli script)
* Files to run snakemake, which include:
  * `workflow` directory, which holds various Snakefile pipelines
  * `snakemake` directory, which defines the pipeline steps and configuration, and shared code between pipelines
