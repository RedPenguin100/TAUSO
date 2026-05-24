"""Binary RNase H1 motif features on the DNA gap of a gapmer ASO.

These sit alongside :func:`seq_features.toxic_motif_count` — both are simple
"is the motif here" presence checks on an ASO sequence — but the RNase H1
motifs are restricted to the *cleavable* DNA gap (the longest consecutive
``d`` stretch in the chemical pattern). Sugar-modified flanks are excluded by
design: RNase H1 does not engage 2'-MOE / cEt / LNA positions, so a TCCC
sitting in the wing is not a cleavage motif.

Strand convention
-----------------
Motifs are evaluated on the **DNA ASO gap** sequence (5'→3'), to match how
Kiełpiński et al. 2017 (PMC5728404) present the motifs textually. Note that
the PSSM scoring paths in ``calculate_rnase`` score the **RNA target strand**
(after ``get_antisense``). PSSM features and motif features therefore look at
opposite strands of the same physical site:

  TCCC on the ASO  ≡  GGGA on the cleaved RNA strand.

If you correlate these features, remember the strand frame they're each in.

Motif selection
---------------
Potency (enhance RNase H1 cleavage) — pyrimidine-rich triplets enriched in
the most-cleaved Kiełpiński substrates:

  - ``TCCC``, ``TTCC``, ``TCTC``

Inefficacy (block cleavage) — a homopolymer G run that can fold into a
G-quadruplex, occluding the substrate from the enzyme:

  - ``GGGG``

The list intentionally excludes pure toxicity flags such as ``TGC`` — those
belong with the toxicity features, not the cleavage-efficacy features.

Gap selection
-------------
Only the **longest** consecutive DNA stretch is considered. Mixmers with two
sizeable gaps (e.g. two 7-nt runs) are deliberately out of scope here —
extending the convention to a "second gap" feature opens questions about
weighting, RNase H1 substrate competition, and interaction with the PSSM
windowing that we are not modelling yet.
"""

from __future__ import annotations

import pandas as pd

from ...common.modifications import get_longest_dna_gap
from ...data.consts import CHEMICAL_PATTERN, SEQUENCE

_RNASEH1_MOTIFS: dict[str, str] = {
    "Potency_TCCC": "TCCC",
    "Potency_TTCC": "TTCC",
    "Potency_TCTC": "TCTC",
    "Inefficacy_GGGG": "GGGG",
}


def _motif_in_gap(gap_seq: str, motif: str) -> float:
    return 1.0 if motif in gap_seq else 0.0


def _extract_dna_gap(row) -> str:
    seq = row.get(SEQUENCE)
    chem = row.get(CHEMICAL_PATTERN)
    if not isinstance(seq, str) or not isinstance(chem, str):
        return ""
    start, end, length = get_longest_dna_gap(chem, marker="d")
    if length < 1:
        return ""
    return seq[start:end].upper()


def add_rnaseh1_motif_features(df: pd.DataFrame, out_prefix: str = "RNaseH1_") -> tuple[pd.DataFrame, list[str]]:
    """Binary presence of curated RNase H1 motifs in the longest DNA gap.

    Adds one float column per motif (1.0 if present in the DNA gap, else 0.0)
    and returns ``(df, feature_cols)``. Rows whose chemistry has no DNA stretch
    or whose sequence/chemistry isn't a string get 0.0 for every motif column.
    """
    gap_sequences = df.apply(_extract_dna_gap, axis=1)

    feature_cols: list[str] = []
    for suffix, motif_seq in _RNASEH1_MOTIFS.items():
        col_name = f"{out_prefix}{suffix}"
        feature_cols.append(col_name)
        df[col_name] = gap_sequences.apply(lambda x, m=motif_seq: _motif_in_gap(x, m))

    return df, feature_cols
