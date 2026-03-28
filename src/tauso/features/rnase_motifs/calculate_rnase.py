import pandas as pd

from ...common.modifications import get_longest_dna_gap
from ...data.consts import CHEMICAL_PATTERN, SEQUENCE
from ...util import get_antisense
from .rnase_helpers import (
    compute_rnaseh1_dinucleotide_score,
    rnaseh1_dict,
    scan_constrained_window,
)


def _apply_rnaseh1_dinuc_scoring(
        df: pd.DataFrame, experiments: tuple[str, ...], out_prefix: str
) -> tuple[pd.DataFrame, list[str]]:
    """Core logic for applying RNase H1 dinucleotide scoring with pre-computation."""
    window_size = 9
    feature_cols: list[str] = []

    # Extract sequences and chemical patterns once as lists
    seqs = df[SEQUENCE].tolist()
    chems = df[CHEMICAL_PATTERN].tolist()

    # --- Pre-compute the gap and target RNA once per row ---
    precomputed_data = []
    for seq, chem in zip(seqs, chems):
        if not isinstance(seq, str) or not isinstance(chem, str):
            precomputed_data.append((None, None))
            continue

        # 1. Find Longest Gap
        start, end, length = get_longest_dna_gap(chem, marker="d")

        # 2. Strict Length Cutoff (Must be >= 8)
        if length < 8:
            precomputed_data.append((None, None))
            continue

        # 3. Extract ASO gap and convert to Target RNA perspective
        aso_gap = seq[start:end]
        target_rna_mimic = get_antisense(aso_gap)

        # 4. Valid starts
        valid_starts = range(len(target_rna_mimic) - window_size + 1)
        if not valid_starts:
            precomputed_data.append((None, None))
        else:
            precomputed_data.append((target_rna_mimic, valid_starts))

    # -------------------------------------------------------------

    # Define the row scorer here. Safely takes 'weights' as an argument.
    def score_row(target_rna_mimic, valid_starts, weights) -> float:
        if target_rna_mimic is None:
            return 0.0

        return max(
            compute_rnaseh1_dinucleotide_score(
                target_rna_mimic, weights, window_start=i
            )
            for i in valid_starts
        )

    for exp in experiments:
        weights = rnaseh1_dict(exp)
        col_name = f"{out_prefix}{exp}_dynamic"
        feature_cols.append(col_name)

        # Loop over our pre-computed list, passing weights directly
        df[col_name] = [
            score_row(t_rna, v_starts, weights) for t_rna, v_starts in precomputed_data
        ]

    return df, feature_cols


def add_rnaseh1_scores_dinuc(
        df: pd.DataFrame,
        experiments: tuple[str, ...] = ("R4a_dinuc", "R4b_dinuc", "R7_dinuc"),
        out_prefix: str = "RNaseH1_score_dinucleotide_",
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
        out_prefix: str = "RNaseH1_Krel_dinucleotide_score_",
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
    seqs = df[SEQUENCE].tolist()
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
        out_prefix: str = "RNaseH1_Krel_score_",
) -> tuple[pd.DataFrame, list[str]]:
    """Calculates RNase H1 SINGLE-NUCLEOTIDE score (LogFC) using best-effort overlap logic."""
    return _apply_rnaseh1_scoring(df, experiments, out_prefix)


def add_rnaseh1_scores_nt(
        df: pd.DataFrame,
        experiments: tuple[str, ...] = ("R4a", "R4b", "R7"),
        out_prefix: str = "RNaseH1_score_",
) -> tuple[pd.DataFrame, list[str]]:
    """Calculates RNase H1 SINGLE-NUCLEOTIDE score (LogFC) using best-effort overlap logic."""
    return _apply_rnaseh1_scoring(df, experiments, out_prefix)


def check_motif_presence(aso_sequence: str, motif: str) -> float:
    """Checks if a motif exists in the sequence (Case Insensitive)."""
    return 1.0 if motif.upper() in aso_sequence.upper() else 0.0


def add_rnaseh1_motif_features(
        df: pd.DataFrame, out_prefix: str = "RNaseH1_"
) -> tuple[pd.DataFrame, list[str]]:
    """
    Extracts the longest DNA gap and checks for specific motifs known to
    enhance (Py-rich) or block (G-quadruplex) RNase H1 efficiency.

    Excludes pure toxicity markers (like TGC).
    """
    # Define motifs: Potency (High Efficacy) and Structural Blocks (Inefficiency)
    rnase_h_motifs = {
        # POSITIVE (High Potency) - Validated in Kiełpiński et al.
        "Potency_TCCC": "TCCC",
        "Potency_TTCC": "TTCC",
        "Potency_TCTC": "TCTC",
        # NEGATIVE - inefficient
        "Inefficacy_GGGG": "GGGG",
    }

    feature_cols = []

    # Pre-calculate the gap strings to avoid re-slicing for every motif
    # We store this temporarily in a Series
    def extract_gap_seq(row):
        seq = row.get(SEQUENCE)
        chem = row.get(CHEMICAL_PATTERN)

        if not isinstance(seq, str) or not isinstance(chem, str):
            return ""

        # Use the robust function to get indices
        start, end, length = get_longest_dna_gap(chem, marker="d")

        if length < 1:
            return ""

        return seq[start:end]

    # Calculate gap sequences once
    gap_sequences = df.apply(extract_gap_seq, axis=1)

    # Generate motif features
    for feature_suffix, motif_seq in rnase_h_motifs.items():
        col_name = f"{out_prefix}{feature_suffix}"
        feature_cols.append(col_name)

        # Apply motif check to the gap sequence series
        df[col_name] = gap_sequences.apply(
            lambda x, m=motif_seq: check_motif_presence(x, m)
        )

    return df, feature_cols
