from collections.abc import Callable, Mapping, Sequence

import pandas as pd

from ...common.modifications import get_longest_dna_gap
from ...data.consts import ASO_SEQUENCE, CHEMICAL_PATTERN
from ...util import get_antisense
from .rnase_helpers import (
    rnaseh1_dict,
    scan_constrained_window,
    scan_constrained_window_dinuc,
)

# The RNase H1 gap window (target RNA, gap start, gap end) for one ASO, or None when
# the ASO has no scorable DNA gap (missing sequence/chemistry, or a gap shorter than 8 nt).
GapWindow = tuple[str, int, int] | None

# A window scorer: (target_rna, weights, gap_start, gap_end) -> score.
WindowScorer = Callable[[str, Mapping[str, float], int, int], float]

MIN_SCORABLE_GAP_LENGTH = 8


def _extract_gap_windows(dataframe: pd.DataFrame) -> list[GapWindow]:
    """Resolve each ASO's DNA-gap window once, so every experiment can reuse it.

    The window depends only on the ASO sequence and chemistry, not on the RNase H1
    weight set, so it is parsed here rather than repeatedly inside the per-experiment loop.
    """
    gap_windows: list[GapWindow] = []
    for aso_sequence, chemistry_pattern in zip(dataframe[ASO_SEQUENCE], dataframe[CHEMICAL_PATTERN]):
        if not isinstance(aso_sequence, str) or not isinstance(chemistry_pattern, str):
            gap_windows.append(None)
            continue

        gap_start_aso, gap_end_aso, gap_length = get_longest_dna_gap(chemistry_pattern, marker="d")
        if gap_length < MIN_SCORABLE_GAP_LENGTH:
            gap_windows.append(None)
            continue

        target_rna = get_antisense(aso_sequence)
        aso_length = len(aso_sequence)
        # Gap coordinates are on the ASO; the target RNA runs antiparallel, so flip them.
        gap_start_target = aso_length - gap_end_aso
        gap_end_target = aso_length - gap_start_aso
        gap_windows.append((target_rna, gap_start_target, gap_end_target))

    return gap_windows


def _apply_rnaseh1_scoring(
    dataframe: pd.DataFrame,
    experiments: Sequence[str],
    out_prefix: str,
    window_scorer: WindowScorer,
) -> tuple[pd.DataFrame, list[str]]:
    """Add one RNase H1 score column per experiment, scored with ``window_scorer``.

    The single-nucleotide and dinucleotide variants differ only in ``window_scorer``;
    everything else — gap resolution and the per-experiment loop — is shared.
    """
    gap_windows = _extract_gap_windows(dataframe)
    feature_cols = []

    for experiment in experiments:
        weights = rnaseh1_dict(experiment)
        col_name = f"{out_prefix}{experiment}_dynamic"
        feature_cols.append(col_name)

        dataframe[col_name] = [
            0.0 if window is None else window_scorer(window[0], weights, window[1], window[2]) for window in gap_windows
        ]

    return dataframe, feature_cols


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
    return _apply_rnaseh1_scoring(df, experiments, out_prefix, scan_constrained_window_dinuc)


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
    return _apply_rnaseh1_scoring(df, experiments, out_prefix, scan_constrained_window_dinuc)


def add_rnaseh1_scores_krel_nt(
    df: pd.DataFrame,
    experiments: tuple[str, ...] = ("R4a_krel", "R4b_krel", "R7_krel"),
    out_prefix: str = "rnase_krel_score_",
) -> tuple[pd.DataFrame, list[str]]:
    """Calculates RNase H1 SINGLE-NUCLEOTIDE score (LogFC) using best-effort overlap logic."""
    return _apply_rnaseh1_scoring(df, experiments, out_prefix, scan_constrained_window)


def add_rnaseh1_scores_nt(
    df: pd.DataFrame,
    experiments: tuple[str, ...] = ("R4a", "R4b", "R7"),
    out_prefix: str = "rnase_score_",
) -> tuple[pd.DataFrame, list[str]]:
    """Calculates RNase H1 SINGLE-NUCLEOTIDE score (LogFC) using best-effort overlap logic."""
    return _apply_rnaseh1_scoring(df, experiments, out_prefix, scan_constrained_window)
