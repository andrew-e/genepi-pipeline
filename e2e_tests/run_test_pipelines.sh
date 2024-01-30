set -e
cp R/tests/data/snakemake_inputs/disease_progression.json input.json
snakemake --snakefile workflow/disease_progression.smk --profile snakemake/bc4/ -F

cp R/tests/data/snakemake_inputs/compare_gwases.json input.json
snakemake --snakefile workflow/compare_gwases.smk --profile snakemake/bc4/ -F

rm input.json
