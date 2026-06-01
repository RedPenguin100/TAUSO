import pandas as pd

from ...common.modifications import get_longest_dna_gap
from ...data.consts import ASO_SEQUENCE, CHEMICAL_PATTERN
from ...util import get_antisense
from .rnase_helpers import (
    rnaseh1_dict,
    scan_constrained_window,
    scan_constrained_window_dinuc,
)


def _apply_rnaseh1_dinuc_scoring(
    df: pd.DataFrame, experiments: tuple[str, ...], out_prefix: str
) -> tuple[pd.DataFrame, list[str]]:
    """Core logic for applying RNase H1 dinucleotide scoring across experiments."""
    feature_cols = []
    seqs = df[ASO_SEQUENCE].tolist()
    chems = df[CHEMICAL_PATTERN].tolist()

    def score_row(seq: str, chem: str, weights: dict) -> float:
        if not isinstance(seq, str) or not isinstance(chem, str):
            return 0.0

        start_aso, end_aso, gap_len = get_longest_dna_gap(chem, marker="d")
        if gap_len < 8:
            return 0.0

        target_rna = get_antisense(seq)
        L = len(seq)

        gap_start_rc = L - end_aso
        gap_end_rc = L - start_aso

        return scan_constrained_window_dinuc(target_rna, weights, gap_start_rc, gap_end_rc)

    for exp in experiments:
        weights = rnaseh1_dict(exp)
        col_name = f"{out_prefix}{exp}_dynamic"
        feature_cols.append(col_name)

        df[col_name] = [score_row(s, c, weights) for s, c in zip(seqs, chems)]

    return df, feature_cols


def add_rnaseh1_scores_dinuc(
    df: pd.DataFrame,
    experiments: tuple[str, ...] = ("R4a_dinuc", "R4b_dinuc", "R7_dinuc"),
    out_prefix: str = "rnase_score_dinucleotide_",
) -> tuple[pd.DataFrame, list[str]]:
    """
    Calculates the best RNase H1 DINUCLEOTIDE score strictly within the longest DNA gap.
    Uses the 'logFC' dinucleotide weights (not krel).
    Returns 0.0 if the DNA gap is shorter than 8 nucleotides.
    """
    return _apply_rnaseh1_dinuc_scoring(df, experiments, out_prefix)


def add_rnaseh1_scores_krel_dinuc(
    df: pd.DataFrame,
    experiments: tuple[str, ...] = (
        "R4a_krel_dinuc",
        "R4b_krel_dinuc",
        "R7_krel_dinuc",
    ),
    out_prefix: str = "rnase_krel_dinucleotide_score_",
) -> tuple[pd.DataFrame, list[str]]:
    """
    Calculates the best RNase H1 score strictly within the longest DNA gap.
    Returns 0.0 if the DNA gap is shorter than 8 nucleotides.
    """
    return _apply_rnaseh1_dinuc_scoring(df, experiments, out_prefix)


def _apply_rnaseh1_scoring(
    df: pd.DataFrame, experiments: tuple[str, ...], out_prefix: str
) -> tuple[pd.DataFrame, list[str]]:
    """Core logic for applying RNase H1 scoring across a set of experiments."""
    feature_cols = []
    seqs = df[ASO_SEQUENCE].tolist()
    chems = df[CHEMICAL_PATTERN].tolist()

    # Define the row scorer here to keep the namespace clean.
    # It safely takes 'weights' as an argument, dodging the late-binding trap.
    def score_row(seq: str, chem: str, weights: dict) -> float:
        if not isinstance(seq, str) or not isinstance(chem, str):
            return 0.0

        start_aso, end_aso, gap_len = get_longest_dna_gap(chem, marker="d")
        if gap_len < 8:
            return 0.0

        target_rna = get_antisense(seq)
        L = len(seq)

        gap_start_rc = L - end_aso
        gap_end_rc = L - start_aso

        return scan_constrained_window(target_rna, weights, gap_start_rc, gap_end_rc)

    for exp in experiments:
        weights = rnaseh1_dict(exp)
        col_name = f"{out_prefix}{exp}_dynamic"
        feature_cols.append(col_name)

        # Pass the loop variable 'weights' directly into the function call
        df[col_name] = [score_row(s, c, weights) for s, c in zip(seqs, chems)]

    return df, feature_cols


def add_rnaseh1_scores_krel_nt(
    df: pd.DataFrame,
    experiments: tuple[str, ...] = ("R4a_krel", "R4b_krel", "R7_krel"),
    out_prefix: str = "rnase_krel_score_",
) -> tuple[pd.DataFrame, list[str]]:
    """Calculates RNase H1 SINGLE-NUCLEOTIDE score (LogFC) using best-effort overlap logic."""
    return _apply_rnaseh1_scoring(df, experiments, out_prefix)


def add_rnaseh1_scores_nt(
    df: pd.DataFrame,
    experiments: tuple[str, ...] = ("R4a", "R4b", "R7"),
    out_prefix: str = "rnase_score_",
) -> tuple[pd.DataFrame, list[str]]:
    """Calculates RNase H1 SINGLE-NUCLEOTIDE score (LogFC) using best-effort overlap logic."""
    return _apply_rnaseh1_scoring(df, experiments, out_prefix)
