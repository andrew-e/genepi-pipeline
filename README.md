# genepi-pipeline
Pipeline for genetic epidemiology projects at Univsersity of Bristol


snakemake \
    -j <MAX JOBS> \
    --configfile <YOUR CONFIG> \
    --use-singularity \
    --singularity-args '--bind $HOME' \
    -s Snakefile

