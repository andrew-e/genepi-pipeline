from dotenv import load_dotenv
from types import SimpleNamespace

def format_dir_string(directory):
    if not directory: return None
    return directory + "/" if not directory.endswith('/') else directory

load_dotenv()
user = os.getenv('USER')
slurm_log_directory = f"/user/work/{user}/slurm_logs/"
docker_container = "docker://andrewrrelmore/genepi_pipeline:latest"

default_clump_headers = "CHR F SNP BP P TOTAL NSIG S05 S01 S001 S0001 SP2"
default_columns = dict(SNP="SNP", CHR="CHR", BP="BP", EA="EA", OA="OA", EAF="EAF", P="P", BETA="BETA",
    SE="SE", OR="OR", OR_LB="OR_LB", OR_UB="OR_UB", RSID="RSID", N="N",
    ENSEMBL_ID="ENSEMBL_ID", GENE_NAME="GENE_NAME"
)
default_mandatory_columns = ["CHR", "BP", "EA", "OA"]
effect_options = ["BETA", "OR", "Z"]
build_options = ["GRCh38", "GRCh37", "GRCh36"]

p_options = ["P", "LOG_P"]
beta_and_or_options = [
    ["BETA", "SE"],
    ["OR", "OR_LB", "OR_UB"]
]
default_output_options = SimpleNamespace(**{
    "build": "GRCh37",
    "columns": SimpleNamespace(**default_columns)
})

if not os.getenv("DATA_DIR") or not os.getenv("RESULTS_DIR"):
    raise ValueError("Please populate DATA_DIR and RESULTS_DIR in the .env file provided")
if not os.getenv("RDFS_DIR"):
    print("Please populate RDFS_DIR in .env if you want the generated files to be automatically copied to RDFS")

DATA_DIR = format_dir_string(os.getenv('DATA_DIR'))
RESULTS_DIR = format_dir_string(os.getenv('RESULTS_DIR'))
RDFS_DIR = format_dir_string(os.getenv('RDFS_DIR'))
THOUSAND_GENOMES_DIR = format_dir_string(os.getenv('GENOMIC_DATA_DIR') + "/1000genomes")
LDSC_DIR = format_dir_string(os.getenv('LDSC_DIR'))
QTL_TOP_HITS_DIR = format_dir_string(os.getenv('QTL_TOP_HITS'))

if RDFS_DIR and not RDFS_DIR.endswith("working/"):
    raise ValueError("Please ensure RDFS_DIR ends with working/ to ensure the data gets copied to the correct place")
