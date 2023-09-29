include: "../snakemake/common.smk"
singularity: docker_container

pipeline = read_json_into_object("input.json")
mandatory_gwas_columns = ["CHR", "BP", "BETA", "SE", "P", "EA", "OA", "EAF"]

onstart:
    print("##### Bluepint for Collider Bias Correction Pipeline #####")

for g in pipeline.gwases:
    g.standardised_gwas = standardised_gwas_name(g.file)
    g.qtl_results = RESULTS_DIR + "mr/" + file_prefix(g.file) + "_" + pipeline.dataset + "_comparison.tsv.gz"

rule all:
    input: qtl_results

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



rule run_mr_against_pqtl_datasets:
    input:
        gwases = [g.standardised_gwas for g in pipeline.gwases],
        ancestries = [g.ancestry for g in pipeline.gwases],
    output: [g.qtl_results for g in pipeline.gwases]
    shell:
        """
        Rscript run_mr_against_qtl_data.r --gwas_filenames {input.gwases} \
            --ancestries {input.ancestries} \
            --dataset {pipeline.dataset} \
            --output {output}
        """

