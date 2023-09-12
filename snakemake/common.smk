import os
from datetime import datetime
from dateutil.relativedelta import relativedelta
from pathlib import Path
from dotenv import load_dotenv
import time

load_dotenv()
user = os.getenv('USER')
slurm_log_directory = f"/user/work/{user}/slurm_logs/"
docker_container = "docker://andrewrrelmore/genepi_pipeline:test"

def validate_gwases(files):
    if not all(Path(file).is_file() for file in files):
        raise ValueError(f"Error: one of the specified files does not exist")

def validate_ancestries(ancestries):
    allowed_ancestries = ["EUR", "EAS", "AFR", "AMR", "SAS"]
    if not all(ancestry in allowed_ancestries for ancestry in ancestries):
        raise ValueError(f"Please ensure all ancestires are one of these values:\n {allowed_ancestries}")
    if len(ancestries) != len(set(ancestries)):
        raise ValueError(f"Error: All GWAS must have a different ancestry")


def standardised_gwas_name(gwas_name):
    if gwas_name.endswith("_standardised.tsv.gz"):
        return gwas_name
    else:
        return DATA_DIR + "gwas/" + file_prefix(gwas_name) + "_standardised.tsv.gz"

def cleanup_old_slurm_logs():
    if not os.path.isdir(slurm_log_directory): return

    one_month_ago = datetime.now() - relativedelta(months=1)
    for filename in os.listdir(slurm_log_directory):
        file = os.path.join(slurm_log_directory, filename)
        file_timestamp = datetime.utcfromtimestamp(os.stat(file).st_mtime)
        if file_timestamp < one_month_ago: os.remove(file)

def input_looper(wildcards):
    print(wildcards)
    n = int(wildcards.n)
    if n == 1:
        return "test/%s.txt" % wildcards.sample
    elif n > 1:
        return "loop%d/%s.txt" % (n-1, wildcards.sample)
    else:
        raise ValueError("loop numbers must be 1 or greater: received %s" % wildcards.n)

def file_prefix(filename):
    stem = Path(rf"{filename}").stem
    return stem.split('.')[0]

def format_dir_string(directory):
    if not directory: return None
    return directory + "/" if not directory.endswith('/') else directory


def turn_results_dict_into_rmd_input(results_dict):
    return ' '.join(['%s=%s' % (key, value) for (key, value) in results_dict.items()])


def copy_data_to_rdfs(files_created):
    if RDFS_DIR is not None:
        for file_created in files_created:
            rdfs_file = None
            if file_created.startswith(DATA_DIR):
                rdfs_file = file_created.replace(DATA_DIR, RDFS_DIR + "/data/")
            elif file_created.startswith(RESULTS_DIR):
                rdfs_file = file_created.replace(RESULTS_DIR, RDFS_DIR + "/results/")

            if rdfs_file:
                os.system(f"install -D {file_created} {rdfs_file}")
        print(f"Files successfully copied to {RDFS_DIR}")


def onsuccess(files_created, results_file=None):
    print("\nWorkflow finished, no errors.  List of created files:")
    print(*files_created, sep='\n')

    if results_file:
        print("PLEASE SEE THIS HTML FOR A SUMMARY OF RESULTS: ", results_file)
        print(f"scp {user}@bc4login1.acrc.bris.ac.uk:{results_file} .")

    copy_data_to_rdfs(files_created)


def onerror_message():
    last_log = subprocess.check_output(f"ls -t {slurm_log_directory} | head -n1", shell=True, universal_newlines=True)
    print("\nThere was an error in the pipeline, please check the last written slurm log to see the error:")
    print(slurm_log_directory + last_log)


if not os.getenv("DATA_DIR") or not os.getenv("RESULTS_DIR"):
    raise ValueError("Please populate DATA_DIR and RESULTS_DIR in the .env file provided")
if not os.getenv("RDFS_DIR"):
    print("Please populate RDFS_DIR in .env if you want the generated files to be automatically copied to RDFS")

if not os.getenv("RDFS_DIR"):
    print("Please populate RDFS_DIR in .env if you want the generated files to be automatically copied to RDFS")

DATA_DIR = format_dir_string(os.getenv('DATA_DIR'))
RESULTS_DIR = format_dir_string(os.getenv('RESULTS_DIR'))
RDFS_DIR = format_dir_string(os.getenv('RDFS_DIR'))

if RDFS_DIR and RDFS_DIR.endswith("working/"):
    raise ValueError("Please ensure RDFS_DIR ends with working/ to ensure the data gets copied to the correct place")

cleanup_old_slurm_logs()
