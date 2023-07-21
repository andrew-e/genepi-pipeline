import os
from pathlib import Path
from dotenv import load_dotenv

load_dotenv()

if not os.getenv("DATA_DIR") or not os.getenv("RESULTS_DIR"):
    raise ValueError("Please populate DATA_DIR and RESULTS_DIR in the .env file provided")

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
    return directory + "/" if not directory.endswith('/') else directory

DATA_DIR = format_dir_string(os.getenv('DATA_DIR'))
RESULTS_DIR = format_dir_string(os.getenv('RESULTS_DIR'))
SCRIPTS_DIR = format_dir_string(os.getenv('SCRIPTS_DIR'))
