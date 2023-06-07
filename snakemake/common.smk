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
