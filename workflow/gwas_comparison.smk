#Basically, turn this lab book into a pipeline:
#https://explodecomputer.github.io/lab-book/posts/2023-06-07-cross-group-effect-comparison/

gwases = [
    {
        "gwas": f"/user/work/{user}/test_data/test_data_no_rsid.tsv.gz",
        "ancestry": "EUR"
    },
    {
        "gwas": f"/user/work/{user}/test_data/subsequent.tsv",
        "ancestry": "EUR"
    }
]

include: "../snakemake/common.smk"
singularity: docker_container

onstart:
    print("##### GWAS Comparison Pipeline #####")
    validate_gwases(g['gwas'] for g in gwases)
    validate_ancestries(g['ancestry'] for g in gwases)

expected_vs_observed_results = RESULTS_DIR + "/expected_vs_observed_outcomes.tsv"
expected_vs_observed_variants = RESULTS_DIR + "/expected_vs_observed_variants.tsv"

for g in gwases:
    g["standardised_gwas"] = DATA_DIR + "gwas/" + file_prefix(g["gwas"]) + "_standardised.tsv.gz"
    g["clumped_snp_prefix"] = DATA_DIR + "clumped_snps/" + file_prefix(g["gwas"])
    g["clumped_snps"] = g["clumped_snp_prefix"] + ".clumped"

rule all:
    input: expected_vs_observed_results, expected_vs_observed_variants

rule standardise_gwases:
    threads: 4
    resources:
        mem = "16G"
    input: ",".join([g['gwas'] for g in gwases])
    output: ",".join([g['standardised_gwas'] for g in gwases])
    shell:
        """
        Rscript standardise_gwas.r \
            --input_gwas {input} \
            --output_gwas {output} \
            --populate_rsid
        """

rule find_clumped_snps:
    resources:
        mem = "4G"
    input:
        clump_dir = DATA_DIR + "clumped_snps",
        gwases = ",".join([g['standardised_gwas'] for g in gwases])
    output:
        ",".join([g['clumped_snps'] for g in gwases])
    shell:
        """
        mkdir -p {input.clump_dir}
        plink1.9 --bfile /user/work/wt23152/genome_data/1000genomes/{ancestry} \
            --clump {input.gwas} \
            --clump-p1 0.00000005 \
            --clump-snp-field RSID \
            --out {clumped_incidence_prefix}
            # TODO: choose these values later --clump-p2 0.001 --clump-r2 0.3 
            #plink1.9 --clump --file {input.gwas} --bfile {input.bfile} --out {output} --clump-r2 0.001 --clump-kb 10000 --clump-p1 1 --clump-p2 1
        """


rule compare_observed_vs_expected_gwas:
    input:
        gwases = ",".join([g['standardised_gwas'] for g in gwases]),
        clumped_files = ",".join([g['clumped_snps'] for g in gwases])
    output:
        results = expected_vs_observed_results,
        variants = expected_vs_observed_variants
    #TODO: maybe change this to run, have a loop, and call shell from there?  Helper function for that?
    shell:
        """
        Rscript compare_observed_vs_expected_gwas.r  \
            --gwas_filenames {gwases} \
            --clumped_filenames {clumped_files} \
            --result_output {output.results} \
            --variants_output {output.variants}
        """

# rule heterogeneity_between_ancestries:
#     input:
#         gwas = expand(DATA_DIR + "/gwas_{ancestry}.tsv", ancestry=ANCESTRIES),
#         clumped_files = expand(DATA_DIR + "/gwas_{ancestry}.tsv", ancestry=ANCESTRIES)
#     output:
#         hetergeneity = RESULTS_DIR + "heterogeneity_score.tsv",
#         heterogeneity_plot_per_snp = RESULTS_DIR + "/plots/hetergeneity_plot.png",
#         comparsion_heterogeneity_snps = RESULTS_DIR + "/plots/comparsion_heterogeneity_snps.png"
#     shell:
#         """
#         Rscript heterogeneity_between_ancestries.r --gwases {gwases} --clumps {clumped_files} \
#             --plot_filename {heterogeneity_plot_per_snp} --comparison_plot_filename {comparsion_heterogeneity_snps}
#         """

files_created = [
    expected_vs_observed_results,
    expected_vs_observed_variants
]

onsuccess:
    onsuccess(files_created)

onerror:
    onerror_message()
