import os
from pathlib import Path

PROJECT_ROOT = Path(os.path.dirname(os.path.dirname(__file__)))
NOTEBOOK_PATH = PROJECT_ROOT / 'notebooks'

ORIGINAL_CSV = NOTEBOOK_PATH / 'data' / 'asoptimizer_data_original.csv'
PROCESSED_CSV = NOTEBOOK_PATH / 'data' / 'asoptimizer_data_processed.csv'
UPDATED_CSV = NOTEBOOK_PATH / 'data' / 'data_asoptimizer_updated.csv'

ORIGINAL_OLIGO_CSV = NOTEBOOK_PATH / 'competitors' / 'oligo_ai' / "aso_inhibitions_21_08_25_incl_context_w_flank_100_df.csv.gz"
PROCESSED_OLIGO_CSV_GZ = NOTEBOOK_PATH / 'competitors' / 'oligo_ai' /'aso_inhibitions_with_canonical_gene.csv.gz'
OLIGO_CSV_INDEXED = NOTEBOOK_PATH / 'competitors' / 'oligo_ai' / 'aso_inhibitions_with_canonical_gene_indexed.csv'

SAVED_FEATURES = NOTEBOOK_PATH / 'saved_features'
CACHE_DIR = PROJECT_ROOT / 'cache'
