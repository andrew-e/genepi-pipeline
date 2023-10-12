import gzip
import json
import os
import time
import re
from datetime import datetime
from dateutil.relativedelta import relativedelta
from pathlib import Path
from dotenv import load_dotenv
from types import SimpleNamespace

load_dotenv()
user = os.getenv('USER')
slurm_log_directory = f"/user/work/{user}/slurm_logs/"
docker_container = "docker://andrewrrelmore/genepi_pipeline:test"
default_clump_headers = "CHR F SNP BP P TOTAL NSIG S05 S01 S001 S0001 SP2"

default_columns = dict(SNP="SNP", CHR="CHR", BP="BP", EA="EA", OA="OA", EAF="EAF", P="P", BETA="BETA", SE="SE", OR="OR", OR_SE="OR_SE", OR_LB="OR_LB", OR_UB="OR_UB", RSID="RSID", N="N")
default_mandatory_columns = ["CHR", "BP", "P", "EA", "OA", "EAF"]

beta_and_or_options = [
    ["BETA", "SE"],
    ["OR", "OR_LB", "OR_UB"],
    ["OR", "OR_SE"]
]

def read_json_into_object(json_file):
    if not os.path.isfile(".env"):
        raise ValueError("Error: .env file doesn't exist")
    if not os.path.isfile(json_file):
        raise ValueError(f"Error: {json_file} file doesn't exist")

    with open(json_file) as pipeline_input:
        pipeline = json.load(pipeline_input,object_hook=lambda data: SimpleNamespace(**data))
    return pipeline


def resolve_gwas_columns(gwas_file, column_name_map=None, additional_mandatory_columns=[]):
    if column_name_map is None:
        column_name_map = SimpleNamespace()

    column_name_map = vars(column_name_map)
    column_name_map = default_columns | column_name_map

    all_mandatory_columns = list(set(default_mandatory_columns + additional_mandatory_columns))
    mandatory_column_names_in_gwas = [column_name_map[name] for name in all_mandatory_columns]

    if not Path(gwas_file).is_file():
        raise ValueError(f"Error: {gwas_file} does not exist")

    if gwas_file.endswith(".gz"):
        with gzip.open(gwas_file) as f:
            gwas_headers = f.readline().decode("utf-8").strip()
    else:
        with open(gwas_file) as f:
            gwas_headers = str(f.readline()).strip()
    gwas_headers = re.split('\n|,| |\t', gwas_headers)

    missing = set(mandatory_column_names_in_gwas) - set(gwas_headers)
    if len(missing) > 0:
        raise ValueError(f"Error: {gwas_file} doesn't contain {missing}")

    beta_and_or_check = []
    for beta_and_or_option in beta_and_or_options:
        option_column_names = [column_name_map[name] for name in beta_and_or_option]
        missing = set(option_column_names) - set(gwas_headers)
        beta_and_or_check.append(len(missing) > 0)

    if all(beta_and_or_check):
        raise ValueError(f"""Error: {gwas_file} doesn't contain the correct pairings for BETA or OR.
            The options available are: (BETA, SE), (OR, OR_SE), or (OR, OR_LB, OR_UB)""")

    return turn_dict_into_cli_string(column_name_map)


def validate_ancestries(ancestries):
    allowed_ancestries = ["EUR", "EAS", "AFR", "AMR", "SAS"]
    if not all(ancestry in allowed_ancestries for ancestry in ancestries):
        raise ValueError(f"Please ensure all ancestries are one of these values:\n {allowed_ancestries}")


def standardised_gwas_name(gwas_name):
    if gwas_name.endswith("_std.tsv.gz"):
        return gwas_name
    else:
        return DATA_DIR + "gwas/" + file_prefix(gwas_name) + "_std.tsv.gz"


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
    prefix = stem.split('.')[0]
    prefix = re.sub("_std", "", prefix)
    return prefix


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
        print("\n---------------------")
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

if RDFS_DIR and not RDFS_DIR.endswith("working/"):
    raise ValueError("Please ensure RDFS_DIR ends with working/ to ensure the data gets copied to the correct place")

cleanup_old_slurm_logs()
