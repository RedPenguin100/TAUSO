import pandas as pd

from ...common.modifications import get_longest_dna_gap
from ...data.consts import CHEMICAL_PATTERN, SEQUENCE
from ...util import get_antisense
from .rnase_helpers import (
    rnaseh1_dict,
    scan_constrained_window,
    scan_constrained_window_dinuc,
)


def _apply_rnaseh1(
    df: pd.DataFrame,
    experiments: tuple[str, ...],
    out_prefix: str,
    scanner,
) -> tuple[pd.DataFrame, list[str]]:
    """Find the longest DNA gap per row, then score it for each experiment.

    Rows whose chemistry has no DNA gap of length >= 8, or that lack a valid
    sequence/chemistry, get 0.0. The DNA gap is the longest consecutive run of
    ``d`` markers — second-largest stretches (e.g. mixmer designs) are ignored.

    ``scanner`` is :func:`scan_constrained_window` for single-nt PSSMs or
    :func:`scan_constrained_window_dinuc` for the dinucleotide variants; both
    zero out positions outside the gap, so sugar-modified flanks contribute
    nothing to the score.
    """
    feature_cols: list[str] = []
    seqs = df[SEQUENCE].tolist()
    chems = df[CHEMICAL_PATTERN].tolist()

    precomputed: list[tuple[str, int, int] | None] = []
    for seq, chem in zip(seqs, chems):
        if not isinstance(seq, str) or not isinstance(chem, str):
            precomputed.append(None)
            continue
        start_aso, end_aso, gap_len = get_longest_dna_gap(chem, marker="d")
        if gap_len < 8:
            precomputed.append(None)
            continue
        L = len(seq)
        target_rna = get_antisense(seq)
        precomputed.append((target_rna, L - end_aso, L - start_aso))

    for exp in experiments:
        weights = rnaseh1_dict(exp)
        col_name = f"{out_prefix}{exp}_dynamic"
        feature_cols.append(col_name)
        df[col_name] = [
            scanner(r[0], weights, r[1], r[2]) if r is not None else 0.0 for r in precomputed
        ]

    return df, feature_cols


def add_rnaseh1_scores_nt(
    df: pd.DataFrame,
    experiments: tuple[str, ...] = ("R4a", "R4b", "R7"),
    out_prefix: str = "RNaseH1_score_",
) -> tuple[pd.DataFrame, list[str]]:
    """Single-nt PSSM RNase H1 scores (logFC weights), zeroed outside the DNA gap."""
    return _apply_rnaseh1(df, experiments, out_prefix, scan_constrained_window)


def add_rnaseh1_scores_krel_nt(
    df: pd.DataFrame,
    experiments: tuple[str, ...] = ("R4a_krel", "R4b_krel", "R7_krel"),
    out_prefix: str = "RNaseH1_Krel_score_",
) -> tuple[pd.DataFrame, list[str]]:
    """Single-nt PSSM RNase H1 scores (Krel weights), zeroed outside the DNA gap."""
    return _apply_rnaseh1(df, experiments, out_prefix, scan_constrained_window)


def add_rnaseh1_scores_dinuc(
    df: pd.DataFrame,
    experiments: tuple[str, ...] = ("R4a_dinuc", "R4b_dinuc", "R7_dinuc"),
    out_prefix: str = "RNaseH1_score_dinucleotide_",
) -> tuple[pd.DataFrame, list[str]]:
    """Dinucleotide PSSM RNase H1 scores (logFC weights), zeroed outside the DNA gap."""
    return _apply_rnaseh1(df, experiments, out_prefix, scan_constrained_window_dinuc)


def add_rnaseh1_scores_krel_dinuc(
    df: pd.DataFrame,
    experiments: tuple[str, ...] = (
        "R4a_krel_dinuc",
        "R4b_krel_dinuc",
        "R7_krel_dinuc",
    ),
    out_prefix: str = "RNaseH1_Krel_dinucleotide_score_",
) -> tuple[pd.DataFrame, list[str]]:
    """Dinucleotide PSSM RNase H1 scores (Krel weights), zeroed outside the DNA gap."""
    return _apply_rnaseh1(df, experiments, out_prefix, scan_constrained_window_dinuc)


def check_motif_presence(aso_sequence: str, motif: str) -> float:
    """Checks if a motif exists in the sequence (case insensitive)."""
    return 1.0 if motif.upper() in aso_sequence.upper() else 0.0


def add_rnaseh1_motif_features(df: pd.DataFrame, out_prefix: str = "RNaseH1_") -> tuple[pd.DataFrame, list[str]]:
    """
    Extracts the longest DNA gap and checks for specific motifs known to
    enhance (Py-rich) or block (G-quadruplex) RNase H1 efficiency.

    Excludes pure toxicity markers (like TGC).
    """
    rnase_h_motifs = {
        # POSITIVE (high potency) — validated in Kiełpiński et al.
        "Potency_TCCC": "TCCC",
        "Potency_TTCC": "TTCC",
        "Potency_TCTC": "TCTC",
        # NEGATIVE — inefficient
        "Inefficacy_GGGG": "GGGG",
    }

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

    for feature_suffix, motif_seq in rnase_h_motifs.items():
        col_name = f"{out_prefix}{feature_suffix}"
        feature_cols.append(col_name)
        df[col_name] = gap_sequences.apply(lambda x, m=motif_seq: check_motif_presence(x, m))

    return df, feature_cols
