include: "../snakemake/common.smk"
singularity: docker_container

pipeline = parse_pipeline_input("input.json")

onstart:
    print("##### Pipeline to Calculate Slope and Apply Correction on Collider Bias #####")

clump_dir = DATA_DIR + "clumped_snps/"
if not os.path.isdir(clump_dir):
    os.makedirs(clump_dir)

incident_file_prefix = file_prefix(pipeline.incident.file)
subsequent_file_prefix = file_prefix(pipeline.subsequent.file)

pipeline.incident.columns = resolve_gwas_columns(pipeline.incident.file, pipeline.incident.columns)
pipeline.subsequent.columns = resolve_gwas_columns(pipeline.subsequent.file, pipeline.subsequent.columns)
pipeline.incident.std_file = standardised_gwas_name(pipeline.incident.file)
pipeline.subsequent.std_file = standardised_gwas_name(pipeline.subsequent.file)
std_file_pattern = standardised_gwas_name("{prefix}")

setattr(pipeline, incident_file_prefix, pipeline.incident)
setattr(pipeline, subsequent_file_prefix, pipeline.subsequent)

clumped_incidence_prefix = clump_dir + incident_file_prefix
clumped_incidence = clumped_incidence_prefix + ".clumped"
clumped_subsequent_prefix = clump_dir + subsequent_file_prefix
clumped_subsequent = clumped_subsequent_prefix + ".clumped"

collider_bias_results = RESULTS_DIR + "collider_bias/" + subsequent_file_prefix + "_collider_bias_results.tsv"
harmonised_effects = RESULTS_DIR + "collider_bias/" + subsequent_file_prefix + "_harmonised_effects.tsv.gz"
slopehunter_results = RESULTS_DIR + "collider_bias/" + subsequent_file_prefix + "_slopehunter.tsv.gz"

unadjusted_miami_plot = RESULTS_DIR + "plots/" + subsequent_file_prefix + "_miami_plot.png"
slopehunter_adjusted_miami_plot = RESULTS_DIR + "plots/" + file_prefix(slopehunter_results) + "_miami_plot.png"
results_file = RESULTS_DIR + "collider_bias/result_" + incident_file_prefix + "_" + subsequent_file_prefix + ".html"
expected_vs_observed_results = RESULTS_DIR + "collider_bias/expected_vs_observed_outcomes.tsv"
expected_vs_observed_variants = RESULTS_DIR + "collider_bias/expected_vs_observed_variants.tsv"


rule all:
    input: expand(std_file_pattern, prefix=[incident_file_prefix, subsequent_file_prefix]),
        collider_bias_results, slopehunter_results, harmonised_effects, unadjusted_miami_plot,
        slopehunter_adjusted_miami_plot, expected_vs_observed_results, expected_vs_observed_variants, results_file


rule standardise_gwases:
    threads: 8
    resources:
        mem = "72G" if pipeline.populate_rsid == "FULL" else "16G"
    params:
        input_gwas = lambda wildcards: getattr(pipeline, wildcards.prefix).file,
        input_columns = lambda wildcards: getattr(pipeline, wildcards.prefix).columns,
    output: std_file_pattern
    shell:
        """
        Rscript standardise_gwas.r \
            --input_gwas {params.input_gwas} \
            --input_columns {params.input_columns} \
            --output_gwas {output} \
            --populate_rsid {pipeline.populate_rsid} 
        """


rule clump_incidence_gwas:
    resources:
        mem = "4G"
    input:
        gwas = pipeline.incident.std_file
    output:
        clumped_incidence
    shell:
        """
        plink1.9 --bfile {THOUSAND_GENOMES_DIR}{pipeline.ancestry} \
            --clump {input.gwas} \
            --clump-snp-field RSID \
            {pipeline.plink_clump_arguments} \
            --out {clumped_incidence_prefix} || echo "{default_clump_headers}" > {output}
        """

#TODO: ensure that this doesn't fail if no hits are found at that p-value
rule clump_subsequent_gwas:
    resources:
        mem = "4G"
    input:
        gwas = pipeline.subsequent.std_file
    output:
        clumped_subsequent
    shell:
        """
        plink1.9 --bfile {THOUSAND_GENOMES_DIR}{pipeline.ancestry} \
            --clump {input.gwas} \
            --clump-snp-field RSID \
            --clump-p1 0.00000005 \
            --out {clumped_subsequent_prefix} || echo "{default_clump_headers}" > {output}
        """

rule collider_bias_correction:
    threads: 8
    resources:
        mem = "32G"
    input:
        incidence_gwas = pipeline.incident.std_file,
        subsequent_gwas = pipeline.subsequent.std_file,
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
        first_gwas = pipeline.incident.std_file,
        second_gwas = pipeline.subsequent.std_file
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
        first_gwas = pipeline.incident.std_file,
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

rule compare_observed_vs_expected_gwas:
    resources:
        mem = "16G"
    input:
        gwases = [pipeline.incident.std_file, pipeline.subsequent.std_file],
        clumped_files = [clumped_incidence]
    output:
        results = expected_vs_observed_results,
        variants = expected_vs_observed_variants
    shell:
        """
        Rscript compare_observed_vs_expected_gwas.r  \
            --gwas_filenames {input.gwases} \
            --clumped_filenames {input.clumped_files} \
            --result_output {output.results} \
            --variants_output {output.variants}
        """

files_created = {
    "incidence_gwas": pipeline.incident.std_file,
    "subsequent_gwas": pipeline.subsequent.std_file,
    "clumped_snps": clumped_incidence,
    "clumped_subsequent": clumped_subsequent,
    "collider_bias_results": collider_bias_results,
    "harmonised_gwas": harmonised_effects,
    "slopehunter_results": slopehunter_results,
    "unadjuested_miami_plot": unadjusted_miami_plot,
    "slopehunter_adjusted_miami_plot": slopehunter_adjusted_miami_plot,
    "expected_vs_observed": expected_vs_observed_results
}
results_string = turn_dict_into_cli_string(files_created)

rule create_results_file:
    threads: 4
    resources:
        mem = "8G",
    input: list(files_created.values())
    output: results_file
    shell:
        """
        Rscript create_results_file.r \
            --rmd_file /home/R/markdown/collider_bias.rmd \
            --params {results_string} \
            --output_file {output}
        """

onsuccess:
    onsuccess(list(files_created.values()), results_file)

onerror:
    onerror_message()
