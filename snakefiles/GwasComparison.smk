include: "../snakemake/common.smk"
singularity: "docker://andrewrrelmore/genepi_pipeline:develop"

#Basically, turn this lab book into a pipeline:
#https://explodecomputer.github.io/lab-book/posts/2023-06-07-cross-group-effect-comparison/


onstart:
    print("##### Bluepint for Multi Ancestry Pipeline #####") 

ANCESTRIES = ["EUR", "EAS", "AFR", "AMR", "SAS"]

rule clump:
    input:
        gwas: expand(DATA_DIR + "/gwas_{ancestry}.tsv", ancestry=ANCESTRIES),
        bfile: expand(DATA_DIR + "{ancestry}", ancestry=ANCESTRIES),
    output:
        expand(DATA_DIR + "{{ancestry}}.clumped", ancestry=ANCESTRIES)
    shell:
        """
        plink1.9 --clump --file {input.gwas} --bfile {input.bfile} --out {output} --clump-r2 0.001 --clump-kb 10000 --clump-p1 1 --clump-p2 1
        """

rule compare_observed_vs_expected_gwas:
    input:
        gwases: ",".join([DATA_DIR + "/gwas_" + a + ".tsv" for a in ANCESTRIES])
        clumped_files: ",".join([DATA_DIR + a + ".clumped" for a in ANCESTRIES])
    output:
        RESULTS_DIR + "/expected_vs_observed_outcomes.tsv"
    shell:
        """
        Rscript compare_observed_vs_expected_gwas.r --gwases {gwases} --clumps {clumped_files} --output {output}
        """

rule heterogeneity_between_ancestries:
    input:
        gwas: expand(DATA_DIR + "/gwas_{ancestry}.tsv", ancestry=ANCESTRIES),
        clumped_files: expand(DATA_DIR + "/gwas_{ancestry}.tsv", ancestry=ANCESTRIES)
    output:
        hetergeneity: RESULTS_DIR + "heterogeneity_score.tsv",
        heterogeneity_plot_per_snp: RESULTS_DIR + "/plots/hetergeneity_plot.png"
        comparsion_heterogeneity_snps: RESULTS_DIR + "/plots/comparsion_heterogeneity_snps.png"
    shell:
        """
        Rscript heterogeneity_between_ancestries.r --gwases {gwases} --clumps {clumped_files} --plot_filename {heterogeneity_plot_per_snp} --comparison_plot_filename {comparsion_heterogeneity_snps}
        """

