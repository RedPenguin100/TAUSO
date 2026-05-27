import logging
import multiprocessing

from ..data.consts import *
from ..parallel_utils import make_apply_fn

logger = logging.getLogger(__name__)
from ..features.hybridization.hybridization_features import (
    calculate_cet,
    calculate_dna,
    calculate_lna,
    get_dna_rna_dg,
    get_dna_rna_dg_region,
    get_ps_delta_dg,
    get_ps_dna_rna_dg,
)
from ..features.hybridization.md_weights import get_moe_md_contribution

# Row-wise hybridization features. All dG values are kcal/mol at 37 C, summed 5'->3'.
HYBR_ROWWISE_CALCULATION = {
    # DNA/RNA hybrid baseline and its phosphorothioate-modified counterpart.
    "hybr_dna_rna_dg": lambda row: get_dna_rna_dg(row[SEQUENCE]),
    "hybr_ps_delta_dg": lambda row: get_ps_delta_dg(row[SEQUENCE]),
    "hybr_ps_dna_rna_dg": lambda row: get_ps_dna_rna_dg(row[SEQUENCE]),
    # DNA/DNA duplex (SantaLucia & Hicks 2004).
    "hybr_dna_dna_dg": lambda row: calculate_dna(row[SEQUENCE]),
    # High-affinity sugar deltas.
    "hybr_lna_delta_dg": lambda row: calculate_lna(row[SEQUENCE], row[CHEMICAL_PATTERN]),
    "hybr_cet_delta_dg": lambda row: calculate_cet(row[SEQUENCE], row[CHEMICAL_PATTERN]),
    # 2'-MOE MD contribution over the MOE-bearing dinucleotides (MOE oligos only).
    "hybr_moe_md_gb_dg": lambda row: get_moe_md_contribution(
        row[SEQUENCE], row[CHEMICAL_PATTERN], row[MODIFICATION], simul_type="gb"
    ),
    "hybr_moe_md_pb_dg": lambda row: get_moe_md_contribution(
        row[SEQUENCE], row[CHEMICAL_PATTERN], row[MODIFICATION], simul_type="pb"
    ),
    # DNA/RNA dG split across the gapmer architecture (5' wing / DNA gap / 3' wing).
    "hybr_dna_rna_dg_wing5": lambda row: get_dna_rna_dg_region(row[SEQUENCE], row[CHEMICAL_PATTERN], "wing5"),
    "hybr_dna_rna_dg_gap": lambda row: get_dna_rna_dg_region(row[SEQUENCE], row[CHEMICAL_PATTERN], "gap"),
    "hybr_dna_rna_dg_wing3": lambda row: get_dna_rna_dg_region(row[SEQUENCE], row[CHEMICAL_PATTERN], "wing3"),
}

# Vectorized features derived from the row-wise ones.
HYBR_DERIVED_FEATURES = [
    "hybr_dna_rna_selectivity_dg",  # DNA/DNA minus DNA/RNA: preference for a DNA over an RNA target
    "hybr_dna_rna_dg_per_nt",
    "hybr_dna_dna_dg_per_nt",
    "hybr_ps_dna_rna_dg_per_nt",
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

    if "hybr_dna_rna_selectivity_dg" in features_to_run and _have("hybr_dna_dna_dg", "hybr_dna_rna_dg"):
        all_data["hybr_dna_rna_selectivity_dg"] = all_data["hybr_dna_dna_dg"] - all_data["hybr_dna_rna_dg"]

    seq_len = all_data[SEQUENCE].str.len()
    normalized = {
        "hybr_dna_rna_dg_per_nt": "hybr_dna_rna_dg",
        "hybr_dna_dna_dg_per_nt": "hybr_dna_dna_dg",
        "hybr_ps_dna_rna_dg_per_nt": "hybr_ps_dna_rna_dg",
    }
    for norm_name, source in normalized.items():
        if norm_name in features_to_run and _have(source):
            all_data[norm_name] = all_data[source] / seq_len

    return all_data, features_to_run
