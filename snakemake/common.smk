import os
from dotenv import load_dotenv

load_dotenv()

if not os.getenv("DATA_DIR") or not os.getenv("RESULTS_DIR"):
    raise ValueError("Please populate DATA_DIR and RESULTS_DIR in the .env file provided")


DATA_DIR = os.getenv('DATA_DIR')
RESULTS_DIR= os.getenv('RESULTS_DIR')
SCRIPTS_DIR= os.getenv('SCRIPTS_DIR')
