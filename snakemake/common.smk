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

def parse_pipeline_input(json_file):
    if not os.path.isfile(".env"):
        raise ValueError("Error: .env file doesn't exist")
    if not os.path.isfile(json_file):
        raise ValueError(f"Error: {json_file} file doesn't exist")

    with open(json_file) as pipeline_input:
        pipeline = json.load(pipeline_input,object_hook=lambda data: SimpleNamespace(**data))

    if not hasattr(pipeline, "output"): pipeline.output = default_output_options
    else:
        if not hasattr(pipeline.output, "build"): setattr(pipeline.output, "build", default_output_options.build)
        elif pipeline.output.build not in build_options:
            raise ValueError(f"Error: {pipeline.output.build} is not one of the valid options: {build_options}")

        if not hasattr(pipeline.output, "columns"): setattr(pipeline.output, "columns", default_output_options.columns)
    if not hasattr(pipeline, "populate_rsid"):
        setattr(pipeline, "populate_rsid", False)

    for g in pipeline.gwases:
        if not hasattr(g, "N"): g.N = 0
        if not hasattr(g, "build"): g.build = "GRCh37"
        g.prefix = file_prefix(g.file)
        g.input_columns = resolve_gwas_columns(g.file, g.columns)
        g.output_columns = resolve_gwas_columns(g.file,pipeline.output.columns, check_input_columns=False)
        g.standardised_gwas = standardised_gwas_name(g.file)
        setattr(pipeline,g.prefix,g)
    return pipeline


def resolve_gwas_columns(gwas_file, column_name_map=None, additional_mandatory_columns=[], check_input_columns=True):
    if not bool(column_name_map):
        column_name_map = SimpleNamespace()

    column_name_map = vars(column_name_map)
    column_name_map = default_columns | column_name_map

    all_mandatory_columns = list(set(default_mandatory_columns + additional_mandatory_columns))
    mandatory_column_names_in_gwas = [column_name_map[name] for name in all_mandatory_columns]

    ensure_mandatory_columns_are_present(gwas_file, mandatory_column_names_in_gwas, column_name_map, check_input_columns)

    cli_string = turn_dict_into_cli_string(column_name_map)
    return cli_string


def ensure_mandatory_columns_are_present(gwas_file, mandatory_column_names_in_gwas, column_name_map, check_input_columns):
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

        p_option_names = [column_name_map.get(name) for name in p_options]
        missing = set(p_option_names) - set(gwas_headers)
        if len(missing) == len(p_option_names):
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

def cleanup_old_slurm_logs():
    if not os.path.isdir(slurm_log_directory): return

    one_month_ago = datetime.now() - relativedelta(months=1)
    files = [f for f in os.listdir(slurm_log_directory) if os.path.isfile(f)]
    print("deleting old logs")

    for filename in files: 
        file = os.path.join(slurm_log_directory, filename)
        file_timestamp = datetime.utcfromtimestamp(os.stat(file).st_mtime)
        if file_timestamp < one_month_ago: os.remove(file)


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
