"""PS-backbone features derived from PS_PATTERN (the ``mod_ps_*`` family).

PS_PATTERN encodes one inter-nucleotide BOND per character ('*' = PS, 'd' = PO), so it has
length ``len(CHEMICAL_PATTERN) - 1`` for a regular oligo. Each bond is attributed to the
region of its 5' nucleotide, matching the convention used in ``get_dna_rna_dg_region``.
"""

import re

import numpy as np
import pandas as pd

from ..common.modifications import get_longest_dna_gap, is_gapmer
from ..data.consts import CHEMICAL_PATTERN, PS_PATTERN

BACKBONE_FEATURES = [
    "mod_ps_po_percentage",
    "mod_ps_end_score",
    "mod_ps_max_consecutive_po",
    "mod_ps_wing5_count",
    "mod_ps_gap_count",
    "mod_ps_wing3_count",
    "mod_ps_frac_mod",
    "mod_ps_frac_dna",
]

# Features that additionally need CHEMICAL_PATTERN (gap/wing placement).
_PLACEMENT_FEATURES = ("mod_ps_wing5_count", "mod_ps_gap_count", "mod_ps_wing3_count")
_INTERACTION_FEATURES = ("mod_ps_frac_mod", "mod_ps_frac_dna")


def max_consecutive_po(ps_pattern):
    if not isinstance(ps_pattern, str):
        return 0
    return max((len(chunk) for chunk in ps_pattern.split("*")), default=0)


def ps_end_score(ps_pattern):
    if not isinstance(ps_pattern, str) or not ps_pattern:
        return 0
    return len(re.match(r"^\**", ps_pattern).group()) + len(re.search(r"\**$", ps_pattern).group())


def ps_placement(chemical_pattern, ps_pattern):
    """PS counts in (wing5, gap, wing3), split by the longest deoxy gap; (0, 0, 0) if no gap."""
    if not isinstance(chemical_pattern, str) or not isinstance(ps_pattern, str):
        return 0, 0, 0
    gap_start, gap_end, gap_len = get_longest_dna_gap(chemical_pattern)
    if gap_len == 0:
        return 0, 0, 0
    return (
        ps_pattern[:gap_start].count("*"),
        ps_pattern[gap_start:gap_end].count("*"),
        ps_pattern[gap_end:].count("*"),
    )


def frac_with_ps(chemical_pattern, ps_pattern, marker_chars):
    """Fraction of positions whose sugar is in ``marker_chars`` that carry a PS bond on their 3' side."""
    if not isinstance(chemical_pattern, str) or not isinstance(ps_pattern, str):
        return 0.0
    # Each PS_PATTERN bond is attributed to the region of its 5' nucleotide
    n = min(len(chemical_pattern), len(ps_pattern))
    qualifying, ps = 0, 0
    for i in range(n):
        if chemical_pattern[i] in marker_chars:
            qualifying += 1
            if ps_pattern[i] == "*":
                ps += 1
    return ps / qualifying if qualifying else 0.0


def _require(df, columns):
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required dependencies in dataframe: {missing}")


def populate_backbone_features(df, features=None):
    """Compute the requested ``mod_ps_*`` backbone features into ``df``; returns ``(df, names)``."""
    names = list(features) if features is not None else list(BACKBONE_FEATURES)
    todo = set(names)

    _require(df, [PS_PATTERN])

    if "mod_ps_max_consecutive_po" in todo:
        df["mod_ps_max_consecutive_po"] = df[PS_PATTERN].apply(max_consecutive_po)

    if "mod_ps_end_score" in todo:
        df["mod_ps_end_score"] = df[PS_PATTERN].apply(ps_end_score)

    if "mod_ps_po_percentage" in todo:
        total_po = df[PS_PATTERN].apply(lambda x: x.count("d") if isinstance(x, str) else 0)
        total_ps = df[PS_PATTERN].apply(lambda x: x.count("*") if isinstance(x, str) else 0)
        backbone_length = total_po + total_ps
        df["mod_ps_po_percentage"] = (total_po / backbone_length.replace(0, pd.NA)).fillna(0)

    need_placement = any(c in todo for c in _PLACEMENT_FEATURES)
    need_interaction = any(c in todo for c in _INTERACTION_FEATURES)

    if need_placement or need_interaction:
        _require(df, [CHEMICAL_PATTERN])

    # mod_ps_wing5/3_count and mod_ps_frac_mod are NaN on non-gapmer rows: their value is 0 by
    # construction there, which would leak zero-variance noise into split selection.
    is_gapmer_series = df[CHEMICAL_PATTERN].map(is_gapmer) if need_placement or need_interaction else None

    if need_placement:
        triples = df.apply(lambda r: ps_placement(r[CHEMICAL_PATTERN], r[PS_PATTERN]), axis=1)
        if "mod_ps_gap_count" in todo:
            # Meaningful for all chemistries (= total PS for all-DNA, 0 for all-modified)
            df["mod_ps_gap_count"] = triples.map(lambda t: t[1]).astype(int)
        if "mod_ps_wing5_count" in todo:
            df["mod_ps_wing5_count"] = triples.map(lambda t: float(t[0])).where(is_gapmer_series, np.nan)
        if "mod_ps_wing3_count" in todo:
            df["mod_ps_wing3_count"] = triples.map(lambda t: float(t[2])).where(is_gapmer_series, np.nan)

    if need_interaction:
        if "mod_ps_frac_mod" in todo:
            # Fraction of high-affinity-sugar positions (M/C/L) that carry a PS bond on their 3' side.
            raw = df.apply(lambda r: frac_with_ps(r[CHEMICAL_PATTERN], r[PS_PATTERN], {"M", "C", "L"}), axis=1).astype(
                float
            )
            df["mod_ps_frac_mod"] = raw.where(is_gapmer_series, np.nan)
        if "mod_ps_frac_dna" in todo:
            # Fraction of deoxy ('d') positions that carry a PS bond on their 3' side.
            df["mod_ps_frac_dna"] = df.apply(
                lambda r: frac_with_ps(r[CHEMICAL_PATTERN], r[PS_PATTERN], {"d"}), axis=1
            ).astype(float)

    return df, names
