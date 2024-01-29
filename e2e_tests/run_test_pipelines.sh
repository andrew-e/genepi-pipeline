conda activate genepi-pipeline || echo "Error: Please ensure your conda environment is initialised and 'genepi-pipeline' env is available" && exit 1
module load apps/singularity/3.8.3

cp R/tests/data/input_compare_gwases.json input.json
snakemake --snakefile workflow/compare_gwases.smk --profile snakemake/bc4/

#cp R/tests/data/input_disease_progression.json input.json
#snakemake --snakefile workflow/disease_progression.smk --profile snakemake/bc4/

rm input.json