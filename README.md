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
* `conda env create --file environment.yml`
* `conda activate genepi-pipeline`
* `module load apps/singularity/3.8.3`

#### Populate the DATA_DIR and RESULTS_DIR and RDFS_DIR environment variables in .env file
These should probably be in your *work* or *scratch* space (`/user/work/$(whoami)/...`)

Note: they can't be your rdfs space, as that volume won't mount to into the singularity container

RDFS_DIR is optional, but very useful.  All generated files will be copied automatically.  Please ensure the path
ends in `working/`, which have the `data` and `results` directories underneath (eg. `/projects/MRC-IEU/research/projects/ieu2/px/xxx/working`)

#### Run the test pipeline
`snakemake --profile snakemake/`

## Usage:

Each pipeline has its own Snakefile in `workflow`.  The standard column naming for GWASes are:

If you do not specify the expected format of the GWAS (ie. plink, bolt, metal), all pipelines will assume the following
headers

|           | SNP | CHR | BP  | BETA | SE  | P   | EA  | OA  | EAF | RSID | OR  | OR_LB | OR_UB |
|-----------|-----|-----|-----|------|-----|-----|-----|-----|-----|:-----|-----|-------|-------|
| Mandatory |     | x   | x   | x    | x   | x   | x   | x   |     |      |     |       |       |

### GWAS Standardisation

All pipelines will standardise each GWAS before running the subsequent steps.  The `SNP` field will be recalculated as

`CHR:POS_EA_OA`, where EA > OA, and the subsequent BETA and EAF will be adjusted accordingly


## How it works:

There are 3 main components to the pipeline
1. Snakemake to define the steps to complete for each pipeline.
2. Docker / Singularity container with intalled languages (R and python), packages, os libraries, and code
3. Slurm: each snakemake step spins up a singularity container inside a slurm job.  Each step can specify different slurm requirements.

## Repository Organisation

* The `Dockerfile` and `docker/` directory hold the information for creating the docker image that the pipeline runs
* `r_scripts` code that can be called from the pipelines, inside the container
  * `r_scripts/` holds the scripts that can be easily called by snakemake (`Rscript example.r --input_file example.txt`)
  * `r_scripts/functions/` holds r code that can be called and reused by any step in the pipeline (accessed by a cli script)
* Files to run snakemake, which include:
  * `workflow` directory, which holds various Snakefile pipelines
  * `snakemake` directory, which defines the pipeline steps and configuration, and shared code between pipelines
