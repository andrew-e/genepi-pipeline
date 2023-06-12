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
        "test.r"
    output:
        work_dir + "/132/results/plots/test.png"
    threads: 1
    shell:
        """
        Rscript /home/scripts/manhattan.r --gwas_filename {input} --manhattan_filename {output.manhattan} --qq_filename {output.qq}
        """
