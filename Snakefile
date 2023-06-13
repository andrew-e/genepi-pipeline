include: "snakemake/common.smk"
singularity: "docker://andrewrrelmore/genepi_pipeline:develop"

onstart:
    print("##### TEST #####") 

rule plinky:
    output: "plink.out"
    shell:
        """
        plink1.9 --help > {output}
        """

rule manhattan_plot:
    resources:
        mem='8G',
        cpus_per_task='4',
        time='01:00:00',
    input:
        DATA_DIR + "/gwases/existing_stroke_gwas.tsv"
    output:
        manhattan = RESULTS_DIR + "/plots/test.png",
        qq = RESULTS_DIR + "/plots/test.qq.png"
    threads: 1
    shell:
        """
        Rscript /home/scripts/manhattan.r --gwas_filename {input} --manhattan_filename {output.manhattan} --qq_filename {output.qq}
        """
