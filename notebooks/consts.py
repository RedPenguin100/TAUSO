import os
from pathlib import Path

PROJECT_ROOT = Path(os.path.dirname(os.path.dirname(__file__)))
NOTEBOOK_PATH = PROJECT_ROOT / 'notebooks'

ORIGINAL_CSV = NOTEBOOK_PATH / 'data' / 'asoptimizer_data_original.csv'
PROCESSED_CSV = NOTEBOOK_PATH / 'data' / 'asoptimizer_data_processed.csv'
UPDATED_CSV = NOTEBOOK_PATH / 'data' / 'data_asoptimizer_updated.csv'


SAVED_FEATURES = NOTEBOOK_PATH / 'saved_features'
CACHE_DIR = PROJECT_ROOT / 'cache'

SEQUENCE = 'Sequence'
INHIBITION = 'Inhibition(%)'
CANONICAL_GENE = 'Canonical Gene Name'
CELL_LINE_ORGANISM = 'Cell line organism'
VOLUME = 'ASO_volume(nM)'
CHEMICAL_PATTERN = 'Chemical_Pattern'
TREATMENT_PERIOD = 'Treatment_Period(hours)'
CELL_LINE = 'Cell_line'
TRANSFECTION = 'Transfection'
DENSITY = 'Density(cells/well)'
DENSITY_UPDATED = 'Density(cells_per_well)' # Avoiding /
PRIMER_PROBE = 'Primer_probe_set'
MODIFICATION = 'Modification'
PREMRNA_FOUND = 'pre_mrna_found'
SENSE_START = 'sense_start'
SENSE_LENGTH = 'sense_length'
SENSE_TYPE = 'sense_type'
SENSE_SEQUENCE = 'sense_sequence'
PRE_MRNA_SEQUENCE = 'pre_mrna_sequence'
SENSE_AVG_ACCESSIBILITY = 'sense_avg_accessibility'

SENSE_SEQUENCE = 'sense_sequence'
PRE_MRNA_SEQUENCE = 'pre_mrna_sequence'
