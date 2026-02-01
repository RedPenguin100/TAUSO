import os
from pathlib import Path

PROJECT_ROOT = Path(os.path.dirname(os.path.dirname(__file__)))
NOTEBOOK_PATH = PROJECT_ROOT / 'notebooks'

ORIGINAL_CSV = NOTEBOOK_PATH / 'data' / 'asoptimizer_data_original.csv'
PROCESSED_CSV = NOTEBOOK_PATH / 'data' / 'asoptimizer_data_processed.csv'
UPDATED_CSV = NOTEBOOK_PATH / 'data' / 'data_asoptimizer_updated.csv'


SAVED_FEATURES = NOTEBOOK_PATH / 'saved_features'
CACHE_DIR = PROJECT_ROOT / 'cache'
