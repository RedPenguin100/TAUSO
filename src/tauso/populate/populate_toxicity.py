"""Populates sequence-derived toxicity / liability features (the tox_* family).

Flags distinct ASO safety liabilities (hepatotoxicity, TLR9 immunostimulation,
G-quadruplex aggregation, acute neurotoxicity). Toxicity is a separate axis from efficacy;
whether these correlate with Inhibition(%) is an empirical question. See
tauso.features.sequence.toxicity_features for the per-mechanism citations.
"""

import logging
import time
from typing import Iterable, Optional, Tuple

import pandas as pd

logger = logging.getLogger(__name__)

from ..data.consts import CHEMICAL_PATTERN, SEQUENCE
from ..features.sequence.toxicity_features import (
    cpg_count,
    g4hunter_max,
    gggg_motif_count,
    hepatotox_gap_motif_count,
    max_g_run,
    three_prime_g_fraction,
    three_prime_g_free_length,
    tlr9_human_motif_count,
)
from ..parallel_utils import make_apply_fn
from ..timer import Timer

# Each feature receives (sequence, chemical_pattern); sequence-only ones ignore the pattern.
FEATURE_SPECS: list[tuple[str, callable]] = [
    # Hepatotoxicity: TCC/TGC in the deoxy gap (Burdick et al. 2014).
    ("tox_hepatotox_gap_tcc_tgc", lambda s, p: hepatotox_gap_motif_count(s, p)),
    # TLR9 immunostimulation: CpG dinucleotide count (Krieg et al. 1995).
    ("tox_cpg_count", lambda s, p: cpg_count(s)),
    # TLR9: human-optimal CpG-B hexamer GTCGTT (Hartmann & Krieg 2000); the canonical
    # human-stimulatory motif, also carried by CpG ODN 2006.
    ("tox_tlr9_gtcgtt", lambda s, p: tlr9_human_motif_count(s)),
    # G-quadruplex propensity: max windowed G4Hunter score (Bedrat et al. 2016).
    ("tox_g4hunter_max", lambda s, p: g4hunter_max(s)),
    # Literal 4+ G-tract count — catches isolated GGGG runs that G4Hunter's
    # windowed mean can underweight.
    ("tox_gggg_motif_count", lambda s, p: gggg_motif_count(s)),
    # Longest run of consecutive G's (G-aggregation / quadruplex proxy).
    ("tox_max_g_run", lambda s, p: max_g_run(s)),
    # Acute neurotoxicity: 3'-end G content and 3'-terminal G-free stretch (Hagedorn et al. 2022).
    ("tox_3prime_g_fraction", lambda s, p: three_prime_g_fraction(s)),
    ("tox_3prime_g_free_len", lambda s, p: three_prime_g_free_length(s)),
]


def calc_feature(df: pd.DataFrame, col_name: str, func, cpus: int = 1, verbose=False) -> None:
    """Computes a feature from the sequence and chemical pattern (parallel if cpus > 1)."""
    start_time = time.time()
    logger.debug("Starting %s...", col_name)
    apply_fn = make_apply_fn(df, n_jobs=cpus, progress_bar=verbose, verbose=0, use_memory_fs=False)
    df[col_name] = apply_fn(lambda row: func(row[SEQUENCE], row[CHEMICAL_PATTERN]), axis=1)
    logger.debug("Finished %s | Time: %.2fs", col_name, time.time() - start_time)


def populate_toxicity_features(
    df,
    features: Optional[Iterable[str]] = None,
    cpus: int = 1,
) -> Tuple:
    available = {name: fn for name, fn in FEATURE_SPECS}
    feature_names = list(features) if features is not None else [name for name, _ in FEATURE_SPECS]

    for name in feature_names:
        with Timer(name):
            calc_feature(df, name, available[name], cpus=cpus)

    return df, feature_names
