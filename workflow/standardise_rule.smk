rule standardise_gwases:
    threads: 8
    resources:
        mem = "72G" if pipeline.populate_rsid == True else "16G"
    params:
        input_gwas = lambda wildcards: getattr(pipeline, wildcards.prefix).file,
        vcf_columns = lambda wildcards: ','.join(getattr(pipeline,wildcards.prefix).columns.values()),
        input_columns = lambda wildcards: getattr(pipeline, wildcards.prefix).input_columns,
        output_columns = lambda wildcards: getattr(pipeline,wildcards.prefix).output_columns
    output: std_file_pattern
    shell:
        """
        INPUT_GWAS={params.input_gwas}
        if [[ {params.input_gwas} =~ ".vcf" ]]; then
            INPUT_GWAS=$(echo "{params.input_gwas}" | sed s/\.vcf/\.tsv/)
            ./vcf_to_tsv.sh {params.input_vcf} {params.vcf_columns} $INPUT_GWAS
        fi

        Rscript standardise_gwas.r \
            --input_gwas $INPUT_GWAS \
            --input_columns {params.input_columns} \
            --output_gwas {output} \
            --populate_rsid {pipeline.populate_rsid} \
            --output_build {pipeline.output.build} \
            --output_effect {pipeline.output.effect} \
            --output_columns {pipeline.output.columns}
        """

