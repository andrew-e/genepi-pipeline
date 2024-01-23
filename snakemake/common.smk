import gzip
import json
import os
import time
import re
from datetime import datetime
from dateutil.relativedelta import relativedelta
from pathlib import Path
from types import SimpleNamespace
include: "constants.smk"

def read_json_into_object(json_file):
    if not os.path.isfile(".env"):
        raise ValueError("Error: .env file doesn't exist")
    if not os.path.isfile(json_file):
        raise ValueError(f"Error: {json_file} file doesn't exist")

    with open(json_file) as pipeline_input:
        pipeline = json.load(pipeline_input,object_hook=lambda data: SimpleNamespace(**data))

    if not hasattr(pipeline, "output"): setattr(pipeline, "output", default_output_options)
    else:
        if not hasattr(pipeline.output, "effect"): setattr(pipeline.output, "effect",default_output_options["effect"])
        elif pipeline.output.effect not in effect_options:
            raise ValueError(f"Error: {pipeline.output.effect} is not one of the valid options: {effect_options}")

        if not hasattr(pipeline.output, "build"): setattr(pipeline.output, "build",default_output_options["build"])
        elif pipeline.output.build not in build_options:
            raise ValueError(f"Error: {pipeline.output.build} is not one of the valid options: {build_options}")

        if not hasattr(pipeline.output, "columns"): setattr(pipeline.output, "columns", default_output_options["columns"])
    if not hasattr(pipeline, "rsid_map"):
        setattr(pipeline, "rsid_map", "PARTIAL")
    elif hasattr(pipeline, "rsid_map") and pipeline.rsid_map not in rsid_map_options:
        raise ValueError(f"Error: {pipeline.rsid_map} is not one of the valid options: {rsid_map_options}")

    return pipeline


def resolve_gwas_columns(gwas_file, column_name_map=None, additional_mandatory_columns=[], check_input_columns=True):
    if column_name_map is None:
        column_name_map = SimpleNamespace()

    column_name_map = vars(column_name_map)
    column_name_map = default_columns | column_name_map

    all_mandatory_columns = list(set(default_mandatory_columns + additional_mandatory_columns))
    mandatory_column_names_in_gwas = [column_name_map[name] for name in all_mandatory_columns]

    ensure_mandatory_columns_are_present(gwas_file, mandatory_column_names_in_gwas, check_input_columns)

    cli_string = turn_dict_into_cli_string(column_name_map)
    return cli_string


def ensure_mandatory_columns_are_present(gwas_file, mandatory_column_names_in_gwas, check_input_columns):
    if not Path(gwas_file).is_file():
        raise ValueError(f"Error: {gwas_file} does not exist")

    if not ".vcf" in gwas_file:
        if gwas_file.endswith(".gz"):
            with gzip.open(gwas_file) as f:
                gwas_headers = f.readline().decode("utf-8").strip()
        else:
            with open(gwas_file) as f:
                gwas_headers = str(f.readline()).strip()
        gwas_headers = re.split('\n|,| |\t',gwas_headers)

        missing = set(mandatory_column_names_in_gwas) - set(gwas_headers)
        if len(missing) > 0 and check_input_columns:
            raise ValueError(f"Error: {gwas_file} doesn't contain {missing}")

        p_option_present = any([p in mandatory_column_names_in_gwas for p in p_options])
        if not p_option_present:
            raise ValueError(f"Error: {gwas_file} doesn't contain a map for P or LOG_P.  Include one please")

        beta_and_or_check = []
        for beta_and_or_option in beta_and_or_options:
            option_column_names = [column_name_map[name] for name in beta_and_or_option]
            missing = set(option_column_names) - set(gwas_headers)
            beta_and_or_check.append(len(missing) > 0)

        if all(beta_and_or_check):
            raise ValueError(f"""Error: {gwas_file} doesn't contain the correct pairings for BETA or OR.
                The options available are: (BETA, SE), or (OR, OR_LB, OR_UB)""")


def validate_ancestries(ancestries):
    allowed_ancestries = ["EUR", "EAS", "AFR", "AMR", "SAS"]
    if not all(ancestry in allowed_ancestries for ancestry in ancestries):
        raise ValueError(f"Please ensure all ancestries are one of these values:\n {allowed_ancestries}")


def standardised_gwas_name(gwas_name):
    if gwas_name.endswith("_std.tsv.gz"):
        return gwas_name
    elif gwas_name.startswith("{") and gwas_name.endswith("}"):
        return DATA_DIR + "gwas/" + gwas_name + "_std.tsv.gz"
    else:
        return DATA_DIR + "gwas/" + file_prefix(gwas_name) + "_std.tsv.gz"

def expanded_gwas_list(gwas_files):
    prefixes = [file_prefix(gwas_file) for gwas_file in gwas_files]
    return expand(DATA_DIR + "gwas/{file_prexix}_std.tsv.gz", file_prefix = prefixes)

def cleanup_old_slurm_logs():
    if not os.path.isdir(slurm_log_directory): return

    one_month_ago = datetime.now() - relativedelta(months=1)
    files = [f for f in os.listdir(slurm_log_directory) if os.path.isfile(f)]
    print("deleting old logs")

    for filename in files: 
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


def turn_dict_into_cli_string(results_dict):
    return ','.join(['%s=%s' % (key, value) for (key, value) in results_dict.items()])


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


def onsuccess(files_created=list(), results_file=None):
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


cleanup_old_slurm_logs()
