include: "../snakemake/common.smk"
singularity: docker_container

with open('input.json') as pipeline_input:
    pipeline = json.load(pipeline_input, object_hook=lambda data: SimpleNamespace(**data))

mandatory_gwas_columns = ["CHR", "BP", "BETA", "SE", "P", "EA", "OA", "EAF"]

onstart:
    print("##### GWAS Ancestry Comparison Pipeline #####")


for gwas in pipeline.gwases:
    gwas.columns = resolve_gwas_columns(gwas.file, gwas.columns, mandatory_gwas_columns)
validate_ancestries([g.ancestry for g in pipeline.gwases])


clump_dir = DATA_DIR + "clumped_snps/"
if not os.path.isdir(clump_dir):
    os.makedirs(clump_dir)

for g in pipeline.gwases:
    g.standardised_gwas = standardised_gwas_name(g.file)
    g.clumped_snp_prefix = clump_dir + file_prefix(g.file)
    g.clumped_snps = g.clumped_snp_prefix + ".clumped"

#List of output files
expected_vs_observed_results = RESULTS_DIR + "ancestry_comparison/expected_vs_observed_outcomes.tsv"
expected_vs_observed_variants = RESULTS_DIR + "ancestry_comparison/expected_vs_observed_variants.tsv"
heterogeneity_scores = RESULTS_DIR + "ancestry_comparison/heterogeneity_scores.tsv"
heterogeneity_plot = RESULTS_DIR + "plots/ancestry_heterogeneity_plot.png"
heterogeneity_snp_comparison = RESULTS_DIR + "plots/ancestry_heterogeneity_snp_comparison.png"
results_file = RESULTS_DIR + "ancestry_comparison/result_summary.html"


ancestries = list([g.ancestry for g in pipeline.gwases])
clumped_snp_prefixes = list([g.clumped_snp_prefix for g in pipeline.gwases])

rule all:
    input: expected_vs_observed_results, expected_vs_observed_variants, heterogeneity_scores, heterogeneity_plot, heterogeneity_snp_comparison, results_file

input_columns = ":".join([g.columns for g in pipeline.gwases])
rule standardise_gwases:
    threads: 4
    resources:
        mem = f"{len(pipeline.gwases)*5}G"
    input: [g.file for g in pipeline.gwases]
    output: [g.standardised_gwas for g in pipeline.gwases]
    shell:
        """
        Rscript standardise_gwas.r \
            --input_gwas {input} \
            --output_gwas {output} \
            --input_columns {input_columns} \
            --populate_rsid
        """

rule find_clumped_snps:
    resources:
        mem = "4G"
    input:
        gwases = [g.standardised_gwas for g in pipeline.gwases]
    output: [g.clumped_snps for g in pipeline.gwases]
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
        mem = f"{len(pipeline.gwases)*10}G"
    input:
        gwases = [g.standardised_gwas for g in pipeline.gwases],
        clumped_files = [g.clumped_snps for g in pipeline.gwases]
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
        mem = f"{len(pipeline.gwases)*10}G"
    input:
        gwases = [g.standardised_gwas for g in pipeline.gwases],
        clumped_files = [g.clumped_snps for g in pipeline.gwases]
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
            --rmd_file markdown/ancestry_comparison.rmd \
            --params {results_string} \
            --output_file {output}
        """

onsuccess:
    onsuccess(list(files_created.values()), results_file)

onerror:
    onerror_message()
