include: "../snakemake/common.smk"
singularity: "docker://andrewrrelmore/genepi_pipeline:develop"

output_file = RESULTS_DIR + "/test_output.log"

onstart:
    print("Test pipeline for genepi-pipeline package")

rule all:
    input: output_file

rule is_everything_installed:
    output: temporary(output_file)
    shell:
        """
        plink1.9 --version >> {output}
        """

