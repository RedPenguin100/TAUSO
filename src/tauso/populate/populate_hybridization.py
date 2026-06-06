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
    "hybr_dna_rna_dg": lambda row: get_dna_rna_dg(row[ASO_SEQUENCE]),
    "hybr_ps_delta_dg": lambda row: get_ps_delta_dg(row[ASO_SEQUENCE], row[PS_PATTERN]),
    "hybr_ps_dna_rna_dg": lambda row: get_ps_dna_rna_dg(row[ASO_SEQUENCE], row[PS_PATTERN]),
    # DNA/DNA duplex (SantaLucia & Hicks 2004).
    "hybr_dna_dna_dg": lambda row: calculate_dna(row[ASO_SEQUENCE]),
    # High-affinity sugar deltas.
    "hybr_lna_delta_dg": lambda row: calculate_lna(row[ASO_SEQUENCE], row[CHEMICAL_PATTERN]),
    "hybr_cet_delta_dg": lambda row: calculate_cet(row[ASO_SEQUENCE], row[CHEMICAL_PATTERN]),
    # 2'-MOE MD contribution over the MOE-bearing dinucleotides (MOE oligos only).
    "hybr_moe_md_gb_dg": lambda row: get_moe_md_contribution(
        row[ASO_SEQUENCE], row[CHEMICAL_PATTERN], row[MODIFICATION_STRING], simul_type="gb"
    ),
    "hybr_moe_md_pb_dg": lambda row: get_moe_md_contribution(
        row[ASO_SEQUENCE], row[CHEMICAL_PATTERN], row[MODIFICATION_STRING], simul_type="pb"
    ),
    # DNA/RNA dG split across the gapmer architecture (5' wing / DNA gap / 3' wing).
    "hybr_dna_rna_dg_wing5": lambda row: get_dna_rna_dg_region(row[ASO_SEQUENCE], row[CHEMICAL_PATTERN], "wing5"),
    "hybr_dna_rna_dg_gap": lambda row: get_dna_rna_dg_region(row[ASO_SEQUENCE], row[CHEMICAL_PATTERN], "gap"),
    "hybr_dna_rna_dg_wing3": lambda row: get_dna_rna_dg_region(row[ASO_SEQUENCE], row[CHEMICAL_PATTERN], "wing3"),
}

# Vectorized features derived from the row-wise ones.
HYBR_DERIVED_FEATURES = [
    "hybr_dna_dna_minus_dna_rna_dg",  # DNA/DNA minus DNA/RNA: the sequence-dependent gap between the
    # DNA/DNA reference (native to the LNA/cEt sugar deltas) and the RNA target -- a reference-state
    # bridge that lets the model relate the DNA-referenced sugar weights to the RNA-bound duplex.
    # NOT a DNA-vs-RNA target preference.
    "hybr_dna_rna_dg_per_nt",
    "hybr_dna_dna_dg_per_nt",
    "hybr_ps_dna_rna_dg_per_nt",
    # Wing/gap asymmetry of the DNA/RNA dG partition. 5'/3' end-stability imbalance is
    # implicated in RNase H1 cleavage-site selectivity and end-mediated nuclease protection;
    # the wing-minus-gap pair quantifies how much the high-affinity flanks dominate the
    # central deoxy region thermodynamically. NaN on non-gapmer rows (no wings to compare).
    "hybr_dna_rna_wing_dg_imbalance",  # wing5_dg - wing3_dg (sign = which end is more stable)
    "hybr_dna_rna_wing5_minus_gap_dg",  # wing5_dg - gap_dg
    "hybr_dna_rna_wing3_minus_gap_dg",  # wing3_dg - gap_dg
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

    seq_len = all_data[ASO_SEQUENCE].str.len()
    normalized = {
        "hybr_dna_rna_dg_per_nt": "hybr_dna_rna_dg",
        "hybr_dna_dna_dg_per_nt": "hybr_dna_dna_dg",
        "hybr_ps_dna_rna_dg_per_nt": "hybr_ps_dna_rna_dg",
    }
    for norm_name, source in normalized.items():
        if norm_name in features_to_run and _have(source):
            all_data[norm_name] = all_data[source] / seq_len

    # Wing/gap asymmetry features are NaN on non-gapmer rows (no wings to compare). The
    # routing flag is derived from CHEMICAL_PATTERN so the populate stays standalone.
    import numpy as np

    from ..common.modifications import get_longest_dna_gap

    def _is_gapmer_pat(pat):
        if not isinstance(pat, str) or not pat:
            return False
        start, end, gap_len = get_longest_dna_gap(pat)
        return gap_len > 0 and start > 0 and end < len(pat)

    wing_asym_features = (
        "hybr_dna_rna_wing_dg_imbalance",
        "hybr_dna_rna_wing5_minus_gap_dg",
        "hybr_dna_rna_wing3_minus_gap_dg",
    )
    is_gapmer_series = (
        all_data[CHEMICAL_PATTERN].map(_is_gapmer_pat)
        if any(f in features_to_run for f in wing_asym_features)
        else None
    )

    if "hybr_dna_rna_wing_dg_imbalance" in features_to_run and _have("hybr_dna_rna_dg_wing5", "hybr_dna_rna_dg_wing3"):
        raw = all_data["hybr_dna_rna_dg_wing5"] - all_data["hybr_dna_rna_dg_wing3"]
        all_data["hybr_dna_rna_wing_dg_imbalance"] = raw.where(is_gapmer_series, np.nan)
    if "hybr_dna_rna_wing5_minus_gap_dg" in features_to_run and _have("hybr_dna_rna_dg_wing5", "hybr_dna_rna_dg_gap"):
        raw = all_data["hybr_dna_rna_dg_wing5"] - all_data["hybr_dna_rna_dg_gap"]
        all_data["hybr_dna_rna_wing5_minus_gap_dg"] = raw.where(is_gapmer_series, np.nan)
    if "hybr_dna_rna_wing3_minus_gap_dg" in features_to_run and _have("hybr_dna_rna_dg_wing3", "hybr_dna_rna_dg_gap"):
        raw = all_data["hybr_dna_rna_dg_wing3"] - all_data["hybr_dna_rna_dg_gap"]
        all_data["hybr_dna_rna_wing3_minus_gap_dg"] = raw.where(is_gapmer_series, np.nan)

    return all_data, features_to_run
