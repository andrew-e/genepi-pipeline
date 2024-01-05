include: "../snakemake/common.smk"
singularity: docker_container

pipeline = read_json_into_object("input.json")

onstart:
    print("##### GWAS LD Score Regression and Genetic Correlation Pipeline #####")


for g in pipeline.gwases:
    g.prefix = file_prefix(g.file)
    g.input_columns = resolve_gwas_columns(gwas.file, gwas.columns, ['N'])
    g.standardised_gwas = standardised_gwas_name(g.file)
    setattr(pipeline, g.prefix, g)

ancestries = list([g.ancestry for g in pipeline.gwases])
validate_ancestries(ancestries)
std_file_pattern = standardised_gwas_name("{prefix}")


rule all:
    input: expand(std_file_pattern, [g.prefix for g in pipeline.gwases]),
        expand(std_file_pattern, [g.prefix for a in ancestries])

rule standardise_gwases:
    threads: 8
    resources:
        mem = "72G" if pipeline.rsid_map == "FULL" else "16G"
    params:
        input_gwas = lambda wildcards: getattr(pipeline, wildcards.prefix).file,
        input_columns = lambda wildcards: getattr(pipeline, wildcards.prefix).input_columns,
        output_columns = lambda wildcards: getattr(pipeline,wildcards.prefix).output_columns
    output: std_file_pattern
    shell:
        """
        Rscript standardise_gwas.r \
            --input_gwas {params.input_gwas} \
            --input_columns {params.input_columns} \
            --output_gwas {output} \
            --populate_rsid {pipeline.rsid_map} 
        """


rule calculate_ldsc_and_genetic_correlation:
    resources:
        mem = f"{len(pipeline.gwases)*4}G"
    input:
        gwases = [g.standardised_gwas for g in pipeline.gwases],
    params:
        gwases = lambda wildcards: getattr(pipeline,wildcards.prefix).file,
    output:
        ldsc_result_log = ldsc_results
    shell:
        """
        ./run_ldsc.sh stuff here
        """
onsuccess:
    onsuccess()

onerror:
    onerror_message()
