include: "../snakemake/common.smk"
singularity: docker_container

pipeline = parse_pipeline_input("input.json")

onstart:
    print("##### GWAS LD Score Regression and Genetic Correlation Pipeline #####")


for g in pipeline.gwases:
    g.prefix = file_prefix(g.file)
    g.input_columns = resolve_gwas_columns(g.file, g.columns, ['N'])
    g.standardised_gwas = standardised_gwas_name(g.file)
    g.sumstats = DATA_DIR + "ldsc/" + g.prefix + ".sumstats.gz"
    setattr(pipeline, g.prefix, g)

ancestries = list(set([g.ancestry for g in pipeline.gwases]))
validate_ancestries(ancestries)
ldsc_result_pattern = RESULTS_DIR + "ldsc/results_{ancestry}.log"
std_file_pattern = standardised_gwas_name("{prefix}")

rule all:
    input: expand(std_file_pattern, prefix=[g.prefix for g in pipeline.gwases]), expand(ldsc_result_pattern, ancestry=ancestries)

rule standardise_gwases:
    threads: 8 if pipeline.populate_rsid == "FULL" else 4
    resources:
        mem = "72G" if pipeline.populate_rsid == "FULL" else "64G"
    params:
        input_gwas = lambda wildcards: getattr(pipeline, wildcards.prefix).file,
        input_columns = lambda wildcards: getattr(pipeline, wildcards.prefix).input_columns,
    output: std_file_pattern
    shell:
        """
        Rscript standardise_gwas.r \
            --input_gwas {params.input_gwas} \
            --input_columns {params.input_columns} \
            --output_gwas {output} \
            --populate_rsid {pipeline.populate_rsid} 
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

onsuccess:
    onsuccess()

onerror:
    onerror_message()
