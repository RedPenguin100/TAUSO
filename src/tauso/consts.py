import os
from pathlib import Path

PACKAGE_PATH = Path(os.path.dirname(os.path.dirname(__file__)))

TEST_PATH = PACKAGE_PATH

DATA_PATH = PACKAGE_PATH / 'data'

GFP1_PATH = DATA_PATH / 'gfp1_seq.txt'
GFP_FIRST_EXP_FASTA = DATA_PATH / 'GFP_first_exp.fasta'

# tmp folder - for files that are dumped to disk during calculation
TMP_PATH = PACKAGE_PATH / 'tmp'

# Experiments
EXPERIMENT_RESULTS = PACKAGE_PATH / 'experiment_results'
CACHE_DIR = PACKAGE_PATH / 'cache'

# External
EXTERNAL_PATH = PACKAGE_PATH / 'external'
RISEARCH_PATH = EXTERNAL_PATH / 'risearch'
RISEARCH1_PATH = RISEARCH_PATH / 'RIsearch1'
RISEARCH1_BINARY_PATH = RISEARCH1_PATH  / 'RIsearch'

# General
APP_NAME = "tauso"
