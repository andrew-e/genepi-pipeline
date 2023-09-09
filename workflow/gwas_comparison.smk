include: "../snakemake/common.smk"
singularity: docker_container

########################################
gwases = [
    {
        "gwas": "/user/work/wt23152/test_data/first_ais_eur_mvp.tsv.gz",
        "ancestry": "EUR"
    },
    {
        "gwas": "/user/work/wt23152/test_data/first_ais_his_mvp.tsv.gz",
        "ancestry": "AMR"
    },
    {
        "gwas": "/user/work/wt23152/test_data/first_ais_afr_mvp.tsv.gz",
        "ancestry": "AFR"
    }
]
########################################

onstart:
    print("##### GWAS Comparison Pipeline #####")
    validate_gwases([g['gwas'] for g in gwases])
    validate_ancestries([g['ancestry'] for g in gwases])

#List of output files
expected_vs_observed_results = RESULTS_DIR + "expected_vs_observed_outcomes.tsv"
expected_vs_observed_variants = RESULTS_DIR + "expected_vs_observed_variants.tsv"
heterogeneity_scores = RESULTS_DIR + "heterogeneity_scores.tsv",
heterogeneity_plot = RESULTS_DIR + "plots/heterogeneity_plot.png",
heterogeneity_snp_comparison = RESULTS_DIR + "plots/heterogeneity_snp_comparison.png"

clump_dir = DATA_DIR + "clumped_snps/"
if not os.path.isdir(clump_dir):
   os.makedirs(clump_dir)

for g in gwases:
    g["standardised_gwas"] = DATA_DIR + "gwas/" + file_prefix(g["gwas"]) + "_standardised.tsv.gz"
    g["clumped_snp_prefix"] = clump_dir + file_prefix(g["gwas"])
    g["clumped_snps"] = g["clumped_snp_prefix"] + ".clumped"

ancestries = list([g['ancestry'] for g in gwases])
clumped_snp_prefixes = list([g['clumped_snp_prefix'] for g in gwases])

rule all:
    input: expected_vs_observed_results, expected_vs_observed_variants, heterogeneity_scores, heterogeneity_plot, heterogeneity_snp_comparison

rule standardise_gwases:
    threads: 4
    resources:
        mem = f"{len(gwases)*5}G"
    input: [g['gwas'] for g in gwases]
    output: [g['standardised_gwas'] for g in gwases]
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
        gwases = [g['standardised_gwas'] for g in gwases]
    output: [g['clumped_snps'] for g in gwases]
    shell:
        """
        gwases=({input.gwases})
        ancestries=({ancestries})
        clumped_snp_prefixes=({clumped_snp_prefixes})

        for i in "${{!gwases[@]}}"
        do
            ancestry=${{ancestries[$i]}}
            plink1.9 --bfile /user/work/wt23152/genome_data/1000genomes/$ancestry \
                --clump ${{gwases[$i]}} \
                --clump-p1 0.0000005 \
                --clump-snp-field RSID \
                --out ${{clumped_snp_prefixes[$i]}}
        done
        """


rule compare_observed_vs_expected_gwas:
    resources:
        mem = f"{len(gwases)*10}G"
    input:
        gwases = [g['standardised_gwas'] for g in gwases],
        clumped_files = [g['clumped_snps'] for g in gwases]
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

rule heterogeneity_between_ancestries:
    resources:
        mem = f"{len(gwases)*10}G"
    input:
        gwases = [g['standardised_gwas'] for g in gwases],
        clumped_files = [g['clumped_snps'] for g in gwases]
    output:
        heterogeneity_scores = heterogeneity_scores,
        heterogeneity_plot = heterogeneity_plot,
        heterogeneity_snp_comparison  = heterogeneity_snp_comparison
    shell:
        """
        Rscript calculate_heterogeneity.r \
            --gwas_filenames {input.gwases} \
            --clumped_filenames {input.clumped_files} \
            --ancestry_list {ancestries} \
            --heterogeneity_scores_output {output.heterogeneity_scores} \
            --heterogeneity_plot_output {output.heterogeneity_plot} \
            --heterogeneity_plot_per_snp_output {output.heterogeneity_snp_comparison}
        """

files_created = [
    expected_vs_observed_results,
    expected_vs_observed_variants,
    heterogeneity_scores,
    heterogeneity_plot,
    heterogeneity_snp_comparison
]

onsuccess:
    onsuccess(files_created)

onerror:
    onerror_message()
