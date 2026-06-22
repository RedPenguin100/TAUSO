"""ASOptimizer constants: raw source-column names and data paths.

These columns come straight from the ASOptimizer source CSVs and do not follow
TAUSO's canonical naming, so every column constant is suffixed ``_ASOPT`` to avoid
colliding with the canonical names in ``tauso.data.consts`` (e.g. ``Cell_line`` here
vs ``cell_line`` there). The cleaning pipeline in ``data_cleanup.py`` reads the raw
CSV with these names and writes out canonical TAUSO columns.

Column names verified against ``raw_data/asoptimizer_data_original.csv``.
"""
from notebooks.consts import NOTEBOOK_PATH

# --- data paths ---
ASOPTIMIZER_DATA_PATH = NOTEBOOK_PATH / "data" / "ASOptimizer" / "raw_data"
ORIGINAL_CSV = ASOPTIMIZER_DATA_PATH / "asoptimizer_data_original.csv"
PROCESSED_CSV = ASOPTIMIZER_DATA_PATH / "asoptimizer_data_processed.csv"

# --- raw ASOptimizer source columns ---
ISIS_ASOPT = "ISIS"
TARGET_GENE_ASOPT = "Target_gene"
CELL_LINE_ASOPT = "Cell_line"
DENSITY_ASOPT = "Density(cells/well)"
TRANSFECTION_ASOPT = "Transfection"
VOLUME_ASOPT = "ASO_volume(nM)"
TREATMENT_PERIOD_ASOPT = "Treatment_Period(hours)"
PRIMER_PROBE_ASOPT = "Primer_probe_set"
SEQUENCE_ASOPT = "Sequence"
MODIFICATION_ASOPT = "Modification"
LOCATION_ASOPT = "Location"
CHEMICAL_PATTERN_ASOPT = "Chemical_Pattern"
LINKAGE_ASOPT = "Linkage"
LINKAGE_LOCATION_ASOPT = "Linkage_Location"
SMILES_ASOPT = "Smiles"
INHIBITION_ASOPT = "Inhibition(%)"
