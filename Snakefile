import os
import getpass

scratch_dir = "/user/work/" + getpass.getuser()
account_id = os.system(f"sacctmgr show user withassoc format=account where user={getpass.getuser()} | grep smed")

print(account_id)

singularity_bindings = {
    "data": scratch_dir + "/132/data",
    "results": scratch_dir + "/132/results"
}

singularity_container = "docker://andrewrrelmore/genepi_pipeline:latest"

def run(command):
    singularity_command = "singularity run "
    for key, value in singularity_bindings.items():
        singularity_command += " ".join([" -B", value + ":" + os.getcwd() + "/" + key,""])

    singularity_command += singularity_container + " " + command
    print(singularity_command)

    return singularity_command

onstart:
    print("##### TEST #####\n") 

rule helloSingularity:
    resources:
        mem='8G',
        cpus_per_task='4',
        time='04:00:00',
    input:
        "test.txt"
    output:
        "out.txt"
    threads: 1
    shell:
        run("Rscript {input} > {output}")

