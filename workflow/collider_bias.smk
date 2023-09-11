include: "../snakemake/common.smk"
singularity: docker_container

#FILL IN VARIABLES BELOW
##############################################
ancestry = "EUR"
incidence_gwas = f"/user/work/{user}/test_data/test_data_no_rsid.tsv.gz"
subsequent_gwas = f"/user/work/{user}/test_data/subsequent.tsv"
plink_clumping_arguments = "--clump-p2 0.001 --clump-r2 0.3"
##############################################


onstart:
    print("##### Pipeline to Calculate Slope and Apply Correction on Collider Bias #####")

standardised_incidence_gwas = standardised_gwas_name(incidence_gwas)
standardised_subsequent_gwas = standardised_gwas_name(subsequent_gwas)

clump_dir = DATA_DIR + "clumped_snps/"
if not os.path.isdir(clump_dir):
   os.makedirs(clump_dir)

clumped_incidence_prefix = clump_dir + file_prefix(incidence_gwas)
clumped_incidence = clumped_incidence_prefix + ".clumped"

collider_bias_results = RESULTS_DIR + "collider_bias/" + file_prefix(subsequent_gwas) + "_collider_bias_results.tsv"
harmonised_effects = RESULTS_DIR + "collider_bias/" + file_prefix(subsequent_gwas) + "_harmonised_effects.tsv.gz"
slopehunter_results = RESULTS_DIR + "collider_bias/" + file_prefix(subsequent_gwas) + "_slopehunter.tsv.gz"

unadjusted_miami_plot = RESULTS_DIR + "plots/" + file_prefix(subsequent_gwas) + "_miami_plot.png"
slopehunter_adjusted_miami_plot = RESULTS_DIR + "plots/" + file_prefix(slopehunter_results) + "_miami_plot.png"

files_created = {
    "incidence_gwas": standardised_incidence_gwas,
    "subsequent_gwas": standardised_subsequent_gwas,
    "clumped_snps": clumped_incidence,
    "collider_bias_results": collider_bias_results,
    "harmonised_gwas": harmonised_effects,
    "slopehunter_results": slopehunter_results,
    "unadjuested_miami_plot": unadjusted_miami_plot,
    "slopehunter_adjusted_miami_plot": slopehunter_adjusted_miami_plot
}
results_file = RESULTS_DIR + "collider_bias/result_" + file_prefix(incidence_gwas) + "_" + file_prefix(subsequent_gwas) + ".html"
results_string = turn_results_dict_into_rmd_input(files_created)

rule all:
    input: collider_bias_results, slopehunter_results, harmonised_effects, unadjusted_miami_plot, slopehunter_adjusted_miami_plot

rule standardise_gwases:
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
        gwas = standardised_incidence_gwas
    output:
        clumped_incidence
    shell:
        """
        plink1.9 --bfile /user/work/wt23152/genome_data/1000genomes/{ancestry} \
            --clump {input.gwas} \
            --clump-snp-field RSID \
            {plink_clumping_arguments} \
            --out {clumped_incidence_prefix}
        """

rule collider_bias_correction:
    threads: 4
    resources:
        mem = "24G"
    input:
        incidence_gwas = standardised_incidence_gwas,
        subsequent_gwas = standardised_subsequent_gwas,
        clumped_file = clumped_incidence
    output:
        results = collider_bias_results,
        slophunter_adjusted = slopehunter_results,
        harmonised_effects_results_file = harmonised_effects
    shell:
        """
        Rscript correct_for_collider_bias.r \
            --incidence_gwas {input.incidence_gwas} \
            --subsequent_gwas {input.subsequent_gwas} \
            --clumped_file {input.clumped_file} \
            --collider_bias_results_output {output.results} \
            --harmonised_effects_output {output.harmonised_effects_results_file} \
            --collider_bias_slopehunter_output {output.slophunter_adjusted}
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

rule create_results_file:
    threads: 4
    resources:
        mem = "8G",
    input: files_created.values()
    output: results_file
    shell:
        """
        Rscript create_results_file.r \
            --rmd_file markdown/collider_bias.rmd \
            --params {results_string} \
            --output_file {output}
        """

onsuccess:
    onsuccess(list(files_created), results_file)

onerror:
    onerror_message()
