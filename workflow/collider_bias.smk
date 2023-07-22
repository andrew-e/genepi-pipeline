include: "../snakemake/common.smk"
singularity: "docker://andrewrrelmore/genepi_pipeline:test"

ancestry = "EUR"
incidence_gwas = "/user/work/wt23152/test_data/incidence.tsv"
subsequent_gwas = "/user/work/wt23152/test_data/subsequent.tsv"

onstart:
    print("##### Pipeline to Calculate Slope and Apply Correction on Collider Bias #####")


clumped_incidence = DATA_DIR + "clumped_snps/" + file_prefix(incidence_gwas)
collider_bias_results = RESULTS_DIR + "collier_bias/" + file_prefix(subsequent_gwas) + "_collider_bias_results.tsv"
slopehunter_results= RESULTS_DIR + "collier_bias/" + file_prefix(subsequent_gwas) + "_slopehunter.tsv.gz"
dudbridge_results = RESULTS_DIR + "collier_bias/" + file_prefix(subsequent_gwas) + "_dudbridge.tsv.gz"

rule all:
    input: collider_bias_results, slopehunter_results, dudbridge_results

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
        slophunter_adjusted = slopehunter_results,
        dudbridge_adjusted = dudbridge_results
    shell:
        """
        Rscript /home/r_scripts/correct_for_collider_bias.r \
            --incidence_gwas {input.incidence_gwas} \
            --subsequent_gwas {input.subsequent_gwas} \
            --clumped_file {input.clumped_file} \
            --collider_bias_results_output {output.results} \
            --collider_bias_slopehunter_output {output.slophunter_adjusted} \
            --collider_bias_dudbridge_output {output.dudbridge_adjusted}
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

