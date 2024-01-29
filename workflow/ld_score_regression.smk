include: "../snakemake/common.smk"
singularity: docker_container

pipeline = parse_pipeline_input("input.json")

onstart:
    print("##### GWAS LD Score Regression and Genetic Correlation Pipeline #####")


for g in pipeline.gwases:
    g.sumstats = DATA_DIR + "ldsc/" + g.prefix + ".sumstats.gz"
    setattr(pipeline, g.prefix, g)

ancestries = list(set([g.ancestry for g in pipeline.gwases]))
validate_ancestries(ancestries)
ldsc_result_pattern = RESULTS_DIR + "ldsc/results_{ancestry}.log"
std_file_pattern = standardised_gwas_name("{prefix}")

rule all:
    input: expand(std_file_pattern, prefix=[g.prefix for g in pipeline.gwases]), expand(ldsc_result_pattern, ancestry=ancestries)

include: "standardise_rule.smk"


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

onsuccess:
    onsuccess()

onerror:
    onerror_message()
