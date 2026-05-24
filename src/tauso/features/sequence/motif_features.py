"""Binary RNase H1 motif presence on the DNA gap of a gapmer ASO.

Sister to :func:`seq_features.toxic_motif_count` — both are sequence-pattern
presence checks. Motifs are scanned on the DNA ASO strand (Kiełpiński 2017
convention, PMC5728404). Add new motifs to ``_RNASEH1_MOTIFS`` below.
"""

import pandas as pd

from ...common.modifications import get_longest_dna_gap
from ...data.consts import CHEMICAL_PATTERN, SEQUENCE

# Potency motifs enhance RNase H1 cleavage; inefficacy motifs (G-quadruplex prone) block it.
# Pure toxicity markers (TGC etc.) belong with the toxicity features, not here.
_RNASEH1_MOTIFS = {
    "Potency_TCCC": "TCCC",
    "Potency_TTCC": "TTCC",
    "Potency_TCTC": "TCTC",
    "Inefficacy_GGGG": "GGGG",
}


def check_motif_presence(aso_sequence: str, motif: str) -> float:
    """Checks if a motif exists in the sequence (case insensitive)."""
    return 1.0 if motif.upper() in aso_sequence.upper() else 0.0


def add_rnaseh1_motif_features(df: pd.DataFrame, out_prefix: str = "RNaseH1_") -> tuple[pd.DataFrame, list[str]]:
    """Add binary presence columns for each motif in ``_RNASEH1_MOTIFS``.

    Only the longest DNA stretch is considered. Mixmers with two sizeable
    gaps are out of scope.
    """
    feature_cols = []

    def extract_gap_seq(row):
        seq = row.get(SEQUENCE)
        chem = row.get(CHEMICAL_PATTERN)

        if not isinstance(seq, str) or not isinstance(chem, str):
            return ""

        start, end, length = get_longest_dna_gap(chem, marker="d")

        if length < 1:
            return ""

        return seq[start:end]

    gap_sequences = df.apply(extract_gap_seq, axis=1)

    for feature_suffix, motif_seq in _RNASEH1_MOTIFS.items():
        col_name = f"{out_prefix}{feature_suffix}"
        feature_cols.append(col_name)
        df[col_name] = gap_sequences.apply(lambda x, m=motif_seq: check_motif_presence(x, m))

    return df, feature_cols
