include: "../snakemake/common.smk"
singularity: "docker://andrewrrelmore/genepi_pipeline:test"

ancestry = "EUR"
incidence_gwas = "/user/work/wt23152/test_data/incidence.tsv"
incidence_gwas = "/user/work/wt23152/test_data/test_data.tsv.gz"
subsequent_gwas = "/user/work/wt23152/test_data/subsequent.tsv"

onstart:
    print("##### Bluepint for Collider Bias Correction Pipeline #####") 


clumped_incidence = DATA_DIR + "clumped_snps/" + file_prefix(incidence_gwas)
collider_bias_results = RESULTS_DIR + "collier_bias/" + file_prefix(subsequent_gwas) + "_collider_bias_results.tsv"
adjusted_results = RESULTS_DIR + "collier_bias/" + file_prefix(subsequent_gwas) + "_collider_bias_adjusted.tsv.gz"
#slopehunter_results= RESULTS_DIR + "collier_bias/" + file_prefix(subsequent_gwas) + "_slopehunter.tsv.gz"
#dudbridge_results = RESULTS_DIR + "collier_bias/" + file_prefix(subsequent_gwas) + "_dudbridge.tsv.gz"

rule all:
    input: collider_bias_results, adjusted_results

#rule standardise_gwas:
#    shell:
#        """
#        Rscript /home/r_scripts/standardise_gwas.r --populate_rsid TRUE
#        """

rule clump_incidence_gwas:
    resources:
        mem = 4000
    input:
        clump_dir = DATA_DIR + "clumped_snps",
        gwas = incidence_gwas,
    output:
        clumped_incidence + ".clumped"
    shell:
        """
        mkdir -p {input.clump_dir}
        plink1.9 --bfile /user/work/wt23152/genome_data/1000genomes/{ancestry} --clump {input.gwas} --clump-snp-field RSID --out {clumped_incidence}
            # TODO: choose these values later--clump-p2 0.001 --clump-r2 0.3 
        """

rule collider_bias_correction:
    threads: 4
    resources:
        mem = 16000
    input:
        incidence_gwas = incidence_gwas,
        subsequent_gwas = subsequent_gwas,
        clumped_file = clumped_incidence + ".clumped"
    output:
        results = collider_bias_results,
        adjusted = adjusted_results
    shell:
        """
        Rscript /home/r_scripts/correct_for_collider_bias.r \
            --incidence_gwas {input.incidence_gwas} \
            --subsequent_gwas {input.subsequent_gwas} \
            --clumped_file {input.clumped_file} \
            --collider_bias_results_output {output.results} \
            --collider_bias_adjusted_output {output.adjusted}
        """

#rule compare_and_plot_collider_bias_corrections:
#    input:
#        gwas: expand(DATA_DIR + "gwas_{ancestry}.tsv", ancestry=ANCESTRIES),
#        clumped_files: expand(DATA_DIR + "gwas_{ancestry}.tsv", ancestry=ANCESTRIES)
#    output:
#        miami_plot_or_something: RESULTS_DIR + "plots/comparsion_heterogeneity_snps.png"
#    shell:
#        """
#        Rscript compare_and_plot_collider_bias_corrections.r --gwases {gwases} \
#            --output_stuff {here}
#        """

