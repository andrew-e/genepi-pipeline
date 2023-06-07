include: "snakemake/common.smk"

onstart:
    print("##### TEST #####") 

rule miami_plot:
    resources:
        mem='8G',
        cpus_per_task='4',
        time='04:00:00',
    input:
        "test.r"
    output:
        "out.txt"
    threads: 1
    shell:
        run("Rscript {input} > {output}")

