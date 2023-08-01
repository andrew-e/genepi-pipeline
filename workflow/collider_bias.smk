include: "../snakemake/common.smk"
singularity: "docker://andrewrrelmore/genepi_pipeline:test"

ancestry = "EUR"
incidence_gwas = "/user/work/wt23152/test_data/test_data_no_rsid.tsv.gz"
subsequent_gwas = "/user/work/wt23152/test_data/subsequent.tsv"

onstart:
    print("##### Pipeline to Calculate Slope and Apply Correction on Collider Bias #####")

standardised_incidence_gwas = DATA_DIR + "gwas/" + file_prefix(incidence_gwas) + "_standardised.tsv.gz"
standardised_subsequent_gwas = DATA_DIR + "gwas/" + file_prefix(subsequent_gwas) + "_standardised.tsv.gz"
clumped_incidence = DATA_DIR + "clumped_snps/" + file_prefix(incidence_gwas)

collider_bias_results = RESULTS_DIR + "collier_bias/" + file_prefix(subsequent_gwas) + "_collider_bias_results.tsv"
slopehunter_results= RESULTS_DIR + "collier_bias/" + file_prefix(subsequent_gwas) + "_slopehunter.tsv.gz"
dudbridge_results = RESULTS_DIR + "collier_bias/" + file_prefix(subsequent_gwas) + "_dudbridge.tsv.gz"

unadjusted_miami_plot = RESULTS_DIR + "plots/" + file_prefix(subsequent_gwas) + "_miami_plot.png"
slopehunter_adjusted_miami_plot = RESULTS_DIR + "plots/" + file_prefix(slopehunter_results) + "_miami_plot.png"

files_created = list(standardised_incidence_gwas, standardised_subsequent_gwas, clumped_incidence, collider_bias_results, slopehunter_results, dudbridge_results, unadjusted_miami_plot, slopehunter_adjusted_miami_plot)

rule all:
    input: collider_bias_results, slopehunter_results, dudbridge_results, unadjusted_miami_plot, slopehunter_adjusted_miami_plot 

rule standardise_gwas:
    threads: 4
    resources:
        mem = "16G"
    input:
        incidence = incidence_gwas,
        subsequent = subsequent_gwas
    output:
        incidence = standardised_incidence_gwas,
        subsequent = standardised_subsequent_gwas
    shell:
        """
        Rscript standardise_gwas.r \
            --input_gwas {input.incidence},{input.subsequent} \
            --output_gwas {output.incidence},{output.subsequent} \
            --populate_rsid
        """

rule clump_incidence_gwas:
    resources:
        mem = "4G"
    input:
        clump_dir = DATA_DIR + "clumped_snps",
        gwas = standardised_incidence_gwas
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
        mem = "16G"
    input:
        incidence_gwas = standardised_incidence_gwas,
        subsequent_gwas = standardised_subsequent_gwas,
        clumped_file = clumped_incidence + ".clumped"
    output:
        results = collider_bias_results,
        slophunter_adjusted = slopehunter_results,
        dudbridge_adjusted = dudbridge_results
    shell:
        """
        Rscript correct_for_collider_bias.r \
            --incidence_gwas {input.incidence_gwas} \
            --subsequent_gwas {input.subsequent_gwas} \
            --clumped_file {input.clumped_file} \
            --collider_bias_results_output {output.results} \
            --collider_bias_slopehunter_output {output.slophunter_adjusted} \
            --collider_bias_dudbridge_output {output.dudbridge_adjusted}
        """

rule unadjusted_miami_plot:
    threads: 4
    resources:
        mem = "8G",
        time = "02:00:00"
    input:
        first_gwas = standardised_incidence_gwas,
        second_gwas = standardised_subsequent_gwas
    output: unadjusted_miami_plot
    shell:
        """
        Rscript miami.r \
            --first_gwas {input.first_gwas} \
            --second_gwas {input.second_gwas} \
            --miami_filename {output} \
            --title "Comparing Incidence and Subsequent GWAS"
        """

rule slopehunter_adjusted_miami_plot:
    threads: 4
    resources:
        mem = "8G",
        time = "02:00:00"
    input:
        first_gwas = standardised_incidence_gwas,
        second_gwas = slopehunter_results
    output: slopehunter_adjusted_miami_plot
    shell:
        """
        Rscript miami.r \
            --first_gwas {input.first_gwas} \
            --second_gwas {input.second_gwas} \
            --miami_filename {output} \
            --title "Comparing Incidence and SlopeHunter Adjusted Subsequent GWAS"
        """

#onsuccess:
#    print("Workflow finished, no errors.  List of created files:")
#    print(*files_created, sep='\n')
#    print(log)
