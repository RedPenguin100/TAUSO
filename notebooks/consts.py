import os
from pathlib import Path

PROJECT_ROOT = Path(os.path.dirname(os.path.dirname(__file__)))
NOTEBOOK_PATH = PROJECT_ROOT / "notebooks"

OLIGO_AI_DATA_PATH = NOTEBOOK_PATH / "data" / "OligoAI" / "raw_data"
ORIGINAL_OLIGO_CSV_FLANK_100 = OLIGO_AI_DATA_PATH / "aso_inhibitions_21_08_25_incl_context_w_flank_100_df.csv.gz"
ORIGINAL_OLIGO_CSV = OLIGO_AI_DATA_PATH / "aso_inhibitions_21_08_25_incl_context_w_flank_50_df.csv.withsplit.csv.gz"
ORIGINAL_OLIGO_CSV_WITH_CANONICAL = OLIGO_AI_DATA_PATH / "aso_inhibitions_with_canonical_gene.csv.gz"

# Raw OligoAI data with our indexing so order doesn't matter
OLIGO_CSV_INDEXED = OLIGO_AI_DATA_PATH / "aso_inhibitions_with_canonical_gene_indexed.csv.gz"
# Processing to enrich the data with information and remove faulty data
OLIGO_CSV_PROCESSED = OLIGO_AI_DATA_PATH / "aso_inhibitions_with_canonical_gene_processed.csv.gz"
# Averaging of rows that are exactly the same to avoid bias in the training / testing phase
OLIGO_CSV_PROCESSED_AVERAGED = OLIGO_AI_DATA_PATH / "aso_inhibitions_with_canonical_gene_processed_averaged.csv.gz"

SAVED_FEATURES = NOTEBOOK_PATH / "saved_features"
CACHE_DIR = PROJECT_ROOT / "cache"

PREDICTION_PATH = NOTEBOOK_PATH / 'prediction'