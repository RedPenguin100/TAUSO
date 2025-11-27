import os
from pathlib import Path

PROJECT_ROOT = Path(os.path.dirname(os.path.dirname(__file__)))
NOTEBOOK_PATH = PROJECT_ROOT / 'notebooks'
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
MODIFICATION = 'Modification'
PREMRNA_FOUND = 'pre_mrna_found'
SENSE_START = 'sense_start'
SENSE_LENGTH = 'sense_length'
SENSE_TYPE = 'sense_type'


