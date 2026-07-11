"""Populates sequence-derived toxicity / liability features (the tox_* family).

Flags distinct ASO safety liabilities (TLR9 immunostimulation,
G-quadruplex aggregation, acute neurotoxicity). Toxicity is a separate axis from efficacy;
whether these correlate with Inhibition(%) is an empirical question. See
tauso.features.sequence.toxicity_features for the per-mechanism citations.
"""

import logging
from typing import Iterable, Optional, Tuple

logger = logging.getLogger(__name__)

from ..data.consts import ASO_SEQUENCE, CHEMICAL_PATTERN, PS_PATTERN
from ..features.sequence.toxicity_features import (
    cpg_count,
    cpg_ps_count,
    g4hunter_max,
    gggg_motif_count,
    three_prime_g_fraction,
    three_prime_g_free_length,
    tlr9_human_motif_count,
)
from .feature_runner import run_feature_specs

# Each feature receives (sequence, chemical_pattern, ps_pattern); features that only need a
# subset just ignore the rest. Keeps the dispatcher signature uniform.
FEATURE_SPECS: list[tuple[str, callable]] = [
    # TLR9 immunostimulation: CpG dinucleotide count (Krieg et al. 1995).
    ("tox_cpg_count", lambda s, p, b: cpg_count(s)),
    # TLR9 amplified by PS backbone: CpG dinucleotides whose bond is PS, the modality used
    # clinically (Krieg 2002; Henry et al. 1997). Targets PS-CpG specifically, which is
    # orders of magnitude more immunostimulatory than the same motif on PO.
    ("tox_cpg_ps_count", lambda s, p, b: cpg_ps_count(s, b)),
    # TLR9: human-optimal CpG-B hexamer GTCGTT (Hartmann & Krieg 2000); the canonical
    # human-stimulatory motif, also carried by CpG ODN 2006.
    ("tox_tlr9_gtcgtt", lambda s, p, b: tlr9_human_motif_count(s)),
    # G-quadruplex propensity: max windowed G4Hunter score (Bedrat et al. 2016).
    ("tox_g4hunter_max", lambda s, p, b: g4hunter_max(s)),
    # Literal 4+ G-tract count — catches isolated GGGG runs that G4Hunter's
    # windowed mean can underweight.
    ("tox_gggg_motif_count", lambda s, p, b: gggg_motif_count(s)),
    # Acute neurotoxicity: 3'-end G content and 3'-terminal G-free stretch (Hagedorn et al. 2022).
    ("tox_3prime_g_fraction", lambda s, p, b: three_prime_g_fraction(s)),
    ("tox_3prime_g_free_len", lambda s, p, b: three_prime_g_free_length(s)),
]


def populate_toxicity_features(
    df,
    features: Optional[Iterable[str]] = None,
    cpus: int = 1,
) -> Tuple:
    return run_feature_specs(
        df,
        FEATURE_SPECS,
        df,
        lambda apply_fn, func: apply_fn(
            lambda row: func(row[ASO_SEQUENCE], row[CHEMICAL_PATTERN], row[PS_PATTERN]), axis=1
        ),
        features=features,
        cpus=cpus,
    )
