import logging
import multiprocessing

from ..data.consts import *
from ..parallel_utils import make_apply_fn

logger = logging.getLogger(__name__)
from ..features.hybridization.hybridization_features import (
    calculate_dna,
    get_cet_dna_rna_dg,
    get_cet_wing_dg,
    get_dna_rna_dg,
    get_lna_dna_rna_dg,
    get_lna_wing_dg,
    get_ps_delta_dg,
    get_ps_dna_rna_dg,
)
from ..features.hybridization.md_weights import get_moe_md_contribution

# Row-wise hybridization features. All dG values are kcal/mol at 37 C, summed 5'->3'.
HYBR_ROWWISE_CALCULATION = {
    # DNA/RNA hybrid baseline and its phosphorothioate-modified counterpart.
    "hybr_dna_rna_dg": lambda row: get_dna_rna_dg(row[ASO_SEQUENCE]),
    "hybr_ps_delta_dg": lambda row: get_ps_delta_dg(row[ASO_SEQUENCE], row[PS_PATTERN]),
    "hybr_ps_dna_rna_dg": lambda row: get_ps_dna_rna_dg(row[ASO_SEQUENCE], row[PS_PATTERN]),
    # DNA/DNA duplex (SantaLucia & Hicks 2004).
    "hybr_dna_dna_dg": lambda row: calculate_dna(row[ASO_SEQUENCE]),
    # Whole-oligo modified-vs-RNA affinity per high-affinity sugar (NaN when that chemistry is
    # absent): cEt/LNA = DNA/RNA baseline + sugar increment; 2'-MOE = full MD duplex energy.
    "hybr_cet_dna_rna_dg": lambda row: get_cet_dna_rna_dg(row[ASO_SEQUENCE], row[CHEMICAL_PATTERN]),
    "hybr_lna_dna_rna_dg": lambda row: get_lna_dna_rna_dg(row[ASO_SEQUENCE], row[CHEMICAL_PATTERN]),
    "hybr_moe_md_gb_dg": lambda row: get_moe_md_contribution(
        row[ASO_SEQUENCE], row[CHEMICAL_PATTERN], row[MODIFICATION_STRING], simul_type="gb"
    ),
    "hybr_moe_md_pb_dg": lambda row: get_moe_md_contribution(
        row[ASO_SEQUENCE], row[CHEMICAL_PATTERN], row[MODIFICATION_STRING], simul_type="pb"
    ),
    # Per-wing high-affinity sugar contribution (5' and 3'); NaN when that chemistry is absent.
    # cEt/LNA = nearest-neighbour increment over the wing; 2'-MOE = MD energy over the wing (GB).
    "hybr_cet_wing5_dg": lambda row: get_cet_wing_dg(row[ASO_SEQUENCE], row[CHEMICAL_PATTERN], "wing5"),
    "hybr_cet_wing3_dg": lambda row: get_cet_wing_dg(row[ASO_SEQUENCE], row[CHEMICAL_PATTERN], "wing3"),
    "hybr_lna_wing5_dg": lambda row: get_lna_wing_dg(row[ASO_SEQUENCE], row[CHEMICAL_PATTERN], "wing5"),
    "hybr_lna_wing3_dg": lambda row: get_lna_wing_dg(row[ASO_SEQUENCE], row[CHEMICAL_PATTERN], "wing3"),
    "hybr_moe_wing5_dg": lambda row: get_moe_md_contribution(
        row[ASO_SEQUENCE], row[CHEMICAL_PATTERN], row[MODIFICATION_STRING], simul_type="gb", region="wing5"
    ),
    "hybr_moe_wing3_dg": lambda row: get_moe_md_contribution(
        row[ASO_SEQUENCE], row[CHEMICAL_PATTERN], row[MODIFICATION_STRING], simul_type="gb", region="wing3"
    ),
}

# Vectorized features derived from the row-wise ones.
HYBR_DERIVED_FEATURES = [
    "hybr_dna_dna_minus_dna_rna_dg",  # DNA/DNA minus DNA/RNA: the sequence-dependent gap between the
    # DNA/DNA reference (native to the LNA/cEt sugar deltas) and the RNA target -- a reference-state
    # bridge that lets the model relate the DNA-referenced sugar weights to the RNA-bound duplex.
    # NOT a DNA-vs-RNA target preference.
]

HYBR_FEATURE_TO_CALCULATION = {
    **HYBR_ROWWISE_CALCULATION,
    **{name: None for name in HYBR_DERIVED_FEATURES},
}


def populate_hybridization(df, n_cores=1, features_to_run=None):
    """Populates hybridization features using the HYBR_FEATURE_TO_CALCULATION registry."""
    all_data = df.copy()

    if features_to_run is None:
        features_to_run = list(HYBR_FEATURE_TO_CALCULATION.keys())

    nb_workers = n_cores if n_cores is not None else multiprocessing.cpu_count()
    apply_func = make_apply_fn(all_data, n_jobs=nb_workers)

    # 1. Row-wise calculations (parallelized).
    for feature in features_to_run:
        logic = HYBR_ROWWISE_CALCULATION.get(feature)
        if callable(logic):
            logger.debug("Calculating feature: %s...", feature)
            all_data[feature] = apply_func(logic, axis=1)

    # 2. Vectorized derived features.
    def _have(*cols):
        missing = [c for c in cols if c not in all_data.columns]
        if missing:
            logger.warning("Skipping derived feature; missing dependencies: %s", missing)
            return False
        return True

    if "hybr_dna_dna_minus_dna_rna_dg" in features_to_run and _have("hybr_dna_dna_dg", "hybr_dna_rna_dg"):
        all_data["hybr_dna_dna_minus_dna_rna_dg"] = all_data["hybr_dna_dna_dg"] - all_data["hybr_dna_rna_dg"]

    return all_data, features_to_run
