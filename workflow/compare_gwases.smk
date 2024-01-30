include: "../snakemake/common.smk"
singularity: docker_container

pipeline = parse_pipeline_input("input.json")

onstart:
    print("##### GWAS Comparison Pipeline #####")

clump_dir = DATA_DIR + "clumped_snps/"
if not os.path.isdir(clump_dir):
    os.makedirs(clump_dir)

for g in pipeline.gwases:
    g.clumped_snp_prefix = clump_dir + file_prefix(g.file)
    g.clumped_file = g.clumped_snp_prefix + ".clumped"
    setattr(pipeline, g.prefix, g)

ancestries = list([g.ancestry for g in pipeline.gwases])
validate_ancestries(ancestries)


#List of output files
expected_vs_observed_results = RESULTS_DIR + "gwas_comparison/expected_vs_observed_outcomes.tsv"
expected_vs_observed_variants = RESULTS_DIR + "gwas_comparison/expected_vs_observed_variants.tsv"
heterogeneity_scores = RESULTS_DIR + "gwas_comparison/heterogeneity_scores.tsv"
heterogeneity_plot = RESULTS_DIR + "plots/ancestry_heterogeneity_plot.png"
heterogeneity_snp_comparison = RESULTS_DIR + "plots/ancestry_heterogeneity_snp_comparison.png"
ldsc_results = RESULTS_DIR + "gwas_comparison/ldsc_results.tsv"
results_file = RESULTS_DIR + "gwas_comparison/result_summary.html"


std_file_pattern = standardised_gwas_name("{prefix}")
ldsc_result_pattern = RESULTS_DIR + "ldsc/results_{ancestry}.log"
rule all:
    input: expand(std_file_pattern, prefix=[g.prefix for g in pipeline.gwases]),
           expand(ldsc_result_pattern, ancestry=ancestries),
           expected_vs_observed_results,
           expected_vs_observed_variants,
           heterogeneity_scores,
           heterogeneity_plot,
           heterogeneity_snp_comparison,
           results_file

include: "standardise_rule.smk"

rule find_clumped_snps:
    resources:
        mem = "4G"
    input:
        gwases = [g.standardised_gwas for g in pipeline.gwases]
    output: [g.clumped_file for g in pipeline.gwases]
    params:
        clumped_snp_prefixes = list([g.clumped_snp_prefix for g in pipeline.gwases])
    shell:
        """
        gwases=({input.gwases})
        ancestries=({ancestries})
        clumped_snp_prefixes=({params.clumped_snp_prefixes})

        for i in "${{!gwases[@]}}"
        do
            ancestry=${{ancestries[$i]}}
            plink1.9 --bfile {THOUSAND_GENOMES_DIR}/$ancestry \
                --clump ${{gwases[$i]}} \
                --clump-p1 0.0000005 \
                --clump-snp-field RSID \
                --out ${{clumped_snp_prefixes[$i]}}
        done
        """


rule compare_observed_vs_expected_gwas:
    resources:
        mem = f"{len(pipeline.gwases)*10}G"
    input:
        gwases = [g.standardised_gwas for g in pipeline.gwases],
        clumped_files = [g.clumped_file for g in pipeline.gwases]
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


rule heterogeneity_between_gwases:
    resources:
        mem = f"{len(pipeline.gwases)*10}G"
    input:
        gwases = [g.standardised_gwas for g in pipeline.gwases],
        clumped_files = [g.clumped_file for g in pipeline.gwases]
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


rule calculate_ldsc_and_genetic_correlation:
    resources:
        mem = "8G"
    params:
        gwases = lambda wildcards: ','.join([g.standardised_gwas for g in pipeline.gwases if g.ancestry == wildcards.ancestry]),
        ns = lambda wildcards: ','.join([str(g.N) for g in pipeline.gwases if g.ancestry == wildcards.ancestry]),
        ancestry = lambda wildcards: wildcards.ancestry
    output: ldsc_result_pattern
    shell:
        """
        ./run_ldsc.sh {params.gwases} {params.ns} {params.ancestry} {output}
        """


files_created = {
    "results": expected_vs_observed_results,
    "variants": expected_vs_observed_variants,
    "heterogeneity_scores": heterogeneity_scores,
    "heterogeneity_plot": heterogeneity_plot,
    "heterogeneity_snp_comparison": heterogeneity_snp_comparison
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
            --rmd_file /home/R/markdown/gwas_comparison.rmd \
            --params {results_string} \
            --output_file {output}
        """

onsuccess:
    onsuccess(list(files_created.values()), results_file)

onerror:
    onerror_message()
