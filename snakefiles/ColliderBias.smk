include: "../snakemake/common.smk"
singularity: "docker://andrewrrelmore/genepi_pipeline:develop"

onstart:
    print("##### Bluepint for Collider Bias Correction Pipeline #####") 

INCIDENCE_GWASES = ["incidence.tsv.gz"]
SUBSEQUENT_GWASES = ["subsequent.tsv.gz", "other_subsequent.tsv.gz"]

rule clump_incidence:
    input:
        gwas: expand(DATA_DIR + "{gwas}", gwas=INCIDENCE_GWASES),
        bfile: expand(DATA_DIR + "{ancestry}", ancestry=ANCESTRIES),
    output:
        expand(DATA_DIR + "/clumped_snps/{{gwas}}.clumped", gwas=INCIDENCE_GWASES)
    #TODO: configure this properly for collider bias clumping / pruning
    shell:
        """
        plink1.9 --clump --file {input.gwas} --bfile {input.bfile} --out {output} --clump-r2 0.001 --clump-kb 10000 --clump-p1 1 --clump-p2 1
        """

rule collider_bias_correction:
    input:
        incidence_gwas: expand(DATA_DIR + "{gwas}", gwas=INCIDENCE_GWASES),
        subsequent_gwas: expand(DATA_DIR + "{gwas}", gwas=SUBSEQUENT_GWASES),
        clumped_file: expand(DATA_DIR + "/clumped_snps/{{gwas}}.clumped", gwas=INCIDENCE_GWASES)
    output:
        RESULTS_DIR + "/expected_vs_observed_outcomes.tsv"
    shell:
        """
        Rscript correct_for_collider_bias.r --incidence_gwas {incidence_gwas} --subsequent_gwas {subsequent_gwas} --clumped_file {clumped_file} --output_file {output}
        """

rule compare_and_plot_collider_bias_corrections:
    input:
        gwas: expand(DATA_DIR + "/gwas_{ancestry}.tsv", ancestry=ANCESTRIES),
        clumped_files: expand(DATA_DIR + "/gwas_{ancestry}.tsv", ancestry=ANCESTRIES)
    output:
        miami_plot_or_something: RESULTS_DIR + "/plots/comparsion_heterogeneity_snps.png"
    shell:
        """
        Rscript compare_and_plot_collider_bias_corrections.r --gwases {gwases} --output_stuff {here}
        """

