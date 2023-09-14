import gzip
import json
import os
import time
from datetime import datetime
from dateutil.relativedelta import relativedelta
from pathlib import Path
from dotenv import load_dotenv
from types import SimpleNamespace

load_dotenv()
user = os.getenv('USER')
slurm_log_directory = f"/user/work/{user}/slurm_logs/"
docker_container = "docker://andrewrrelmore/genepi_pipeline:test"
default_columns = dict(SNP="SNP", CHR="CHR", BP="BP", EA="EA", OA="OA", EAF="EAF", P="P", BETA="BETA", SE="SE", OR="OR", OR_LB="OR_LB", OR_UB="OR_UB", RSID="RSID")

def resolve_gwas_columns(gwas_file, columns, mandatory_gwas_columns):
    columns = vars(columns)
    columns = default_columns | columns

    if not Path(gwas_file).is_file():
        raise ValueError(f"Error: {gwas_file} does not exist")

    if gwas_file.endswith(".gz"):
        with gzip.open(gwas_file) as f:
            headers = f.readline().decode("utf-8").strip()
    else:
        with open(gwas_file) as f:
            headers = str(f.readline()).strip()
    headers = re.split('\n|,| |\t', headers)

    for mandatory_column in mandatory_gwas_columns:
        column = columns[mandatory_column]
        if not column in headers:
            raise ValueError(f"Error: {gwas_file} doesn't contain mandatory column {mandatory_column}/{column}")

    return turn_dict_into_cli_string(columns)

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


def turn_dict_into_cli_string(results_dict):
    return ','.join(['%s=%s' % (key, value) for (key, value) in results_dict.items()])

def ensure_mandatory_columns_are_present(gwas, columns=default_columns):
    return 0

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
        print("PLEASE SEE THIS HTML FOR A SUMMARY OF RESULTS:")
        print(f"scp {user}@bc4login1.acrc.bris.ac.uk:{results_file} .")

    copy_data_to_rdfs(files_created)


def onerror_message():
    last_log = subprocess.check_output(f"ls -t {slurm_log_directory} | head -n1", shell=True, universal_newlines=True)
    print("\n---------------------")
    print("There was an error in the pipeline, please check the last written slurm log to see the error:")
    print(slurm_log_directory + last_log)


if not os.getenv("DATA_DIR") or not os.getenv("RESULTS_DIR"):
    raise ValueError("Please populate DATA_DIR and RESULTS_DIR in the .env file provided")
if not os.getenv("RDFS_DIR"):
    print("Please populate RDFS_DIR in .env if you want the generated files to be automatically copied to RDFS")

DATA_DIR = format_dir_string(os.getenv('DATA_DIR'))
RESULTS_DIR = format_dir_string(os.getenv('RESULTS_DIR'))
RDFS_DIR = format_dir_string(os.getenv('RDFS_DIR'))
THOUSAND_GENOMES_DIR = format_dir_string(os.getenv('THOUSAND_GENOMES_DIR'))
LDSC_DIR = format_dir_string(os.getenv('LDSC_DIR'))
QTL_TOP_HITS_DIR = format_dir_string(os.getenv('QTL_TOP_HITS'))


if RDFS_DIR and RDFS_DIR.endswith("working/"):
    raise ValueError("Please ensure RDFS_DIR ends with working/ to ensure the data gets copied to the correct place")

cleanup_old_slurm_logs()
