include: "../snakemake/common.smk"
singularity: "docker://andrewrrelmore/genepi_pipeline:develop"

###SPECIFY INPUTS HERE
gwas = "gwas.tsv.gz"
###SPECIFY INPUTS HERE

onstart:
    print("##### Bluepint for Collider Bias Correction Pipeline #####") 


pqtl_results = RESULTS_DIR + "mr/" + file_prefix(subsequent_gwas) + "_pqtl_comparison.tsv.gz"

rule all:
    input: pqtl_results

rule run_mr_against_pqtl_datasets:
    input: gwas 
    output: pqtl_results
    shell:
        """
        Rscript run_mr_against_pqtl_data.r --gwas {input} --output {output}
        """
