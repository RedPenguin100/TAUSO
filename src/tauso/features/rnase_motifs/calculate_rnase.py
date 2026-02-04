import pandas as pd

from ...data.consts import SEQUENCE, CHEMICAL_PATTERN
from .rnase_helpers import compute_rnaseh1_dinucleotide_score, rnaseh1_dict, \
    get_longest_dna_gap, scan_constrained_window
from ...util import get_antisense

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
    # The Dinucleotide models (R4a_dinuc, etc.) typically use a 9-mer window
    window_size = 9

    feature_cols: list[str] = []

    for exp in experiments:
        weights = rnaseh1_dict(exp)
        col_name = f"{out_prefix}{exp}_dynamic"
        feature_cols.append(col_name)

        def score_row(row) -> float:
            seq = row.get(SEQUENCE)
            chem = row.get(CHEMICAL_PATTERN)

            if not isinstance(seq, str) or not isinstance(chem, str):
                return 0.0

            # 1. Find Longest Gap (ASO coordinates)
            start, end, length = get_longest_dna_gap(chem, marker="d")

            # 2. Strict Length Cutoff (Must be >= 8)
            if length < 8:
                return 0.0

            # 3. Extract ASO gap and convert to Target RNA perspective
            aso_gap = seq[start:end]
            target_rna_mimic = get_antisense(aso_gap)

            # 4. Scan for best score within the TARGET sequence
            valid_starts = range(len(target_rna_mimic) - window_size + 1)

            if not valid_starts:
                return 0.0

            return max(
                compute_rnaseh1_dinucleotide_score(target_rna_mimic, weights, window_start=i)
                for i in valid_starts
            )

        df[col_name] = df.apply(score_row, axis=1)

    return df, feature_cols

def add_rnaseh1_scores_krel_dinuc(
        df: pd.DataFrame,
        experiments: tuple[str, ...] = ("R4a_krel_dinuc", "R4b_krel_dinuc", "R7_krel_dinuc"),
        out_prefix: str = "RNaseH1_Krel_dinucleotide_score_",
) -> tuple[pd.DataFrame, list[str]]:
    """
    Calculates the best RNase H1 score strictly within the longest DNA gap.
    Returns 0.0 if the DNA gap is shorter than 8 nucleotides.
    """
    window_size = 9

    feature_cols: list[str] = []

    for exp in experiments:
        weights = rnaseh1_dict(exp)
        col_name = f"{out_prefix}{exp}_dynamic"
        feature_cols.append(col_name)

        def score_row(row) -> float:
            seq = row.get(SEQUENCE)
            chem = row.get(CHEMICAL_PATTERN)

            if not isinstance(seq, str) or not isinstance(chem, str):
                return 0.0

            # 1. Find Longest Gap
            start, end, length = get_longest_dna_gap(chem, marker="d")

            # 2. Strict Length Cutoff (Must be >= 8)
            if length < 8:
                return 0.0

            # 3. Extract ASO gap and convert to Target RNA perspective
            aso_gap = seq[start:end]
            target_rna_mimic = get_antisense(aso_gap)

            # 4. Scan for best score within the TARGET sequence
            # Note: valid_starts is based on the new target_rna_mimic length (same as gap length)
            valid_starts = range(len(target_rna_mimic) - window_size + 1)

            if not valid_starts:
                return 0.0

            return max(
                compute_rnaseh1_dinucleotide_score(target_rna_mimic, weights, window_start=i)
                for i in valid_starts
            )

        df[col_name] = df.apply(score_row, axis=1)

    return df, feature_cols


def add_rnaseh1_scores_krel_nt(
        df: pd.DataFrame,
        experiments: tuple[str, ...] = ("R4a_krel", "R4b_krel", "R7_krel"),
        out_prefix: str = "RNaseH1_Krel_score_",
) -> tuple[pd.DataFrame, list[str]]:
    """
    Calculates RNase H1 SINGLE-NUCLEOTIDE score using best-effort overlap logic.
    """
    feature_cols: list[str] = []

    for exp in experiments:
        weights = rnaseh1_dict(exp)
        col_name = f"{out_prefix}{exp}_dynamic"
        feature_cols.append(col_name)

        def score_row(row) -> float:
            seq = row.get(SEQUENCE)
            chem = row.get(CHEMICAL_PATTERN)

            if not isinstance(seq, str) or not isinstance(chem, str):
                return 0.0

            # 1. Get DNA Gap indices (ASO coordinates)
            start_aso, end_aso, gap_len = get_longest_dna_gap(chem, marker="d")
            if gap_len < 2: return 0.0

            # 2. Convert to Target RNA perspective
            target_rna = get_antisense(seq)
            L = len(seq)

            # Map ASO gap indices [start, end) to Target indices [rc_start, rc_end)
            # The gap is flipped: End of ASO is Start of Target
            gap_start_rc = L - end_aso
            gap_end_rc = L - start_aso

            # 3. Calculate Score using helper
            return scan_constrained_window(target_rna, weights, gap_start_rc, gap_end_rc)

        df[col_name] = df.apply(score_row, axis=1)

    return df, feature_cols


def add_rnaseh1_scores_nt(
        df: pd.DataFrame,
        experiments: tuple[str, ...] = ("R4a", "R4b", "R7"),
        out_prefix: str = "RNaseH1_score_",
) -> tuple[pd.DataFrame, list[str]]:
    """
    Calculates RNase H1 SINGLE-NUCLEOTIDE score (LogFC) using best-effort overlap logic.
    Replaces hardcoded 'best window' maps with dynamic gap scanning.
    """
    feature_cols: list[str] = []

    for exp in experiments:
        weights = rnaseh1_dict(exp)
        col_name = f"{out_prefix}{exp}_dynamic"
        feature_cols.append(col_name)

        def score_row(row) -> float:
            seq = row.get(SEQUENCE)
            chem = row.get(CHEMICAL_PATTERN)

            if not isinstance(seq, str) or not isinstance(chem, str):
                return 0.0

            # 1. Get DNA Gap indices (ASO coordinates)
            start_aso, end_aso, gap_len = get_longest_dna_gap(chem, marker="d")
            if gap_len < 2: return 0.0

            # 2. Convert to Target RNA perspective
            target_rna = get_antisense(seq)
            L = len(seq)

            # Map ASO gap indices [start, end) to Target indices [rc_start, rc_end)
            # The gap is flipped: End of ASO is Start of Target
            gap_start_rc = L - end_aso
            gap_end_rc = L - start_aso

            # 3. Calculate Score using helper
            # This automatically handles the "R4a/R7" weights using the same
            # logic as the krel function (best effort overlap).
            return scan_constrained_window(target_rna, weights, gap_start_rc, gap_end_rc)

        df[col_name] = df.apply(score_row, axis=1)

    return df, feature_cols


def check_motif_presence(aso_sequence: str, motif: str) -> float:
    """Checks if a motif exists in the sequence (Case Insensitive)."""
    return 1.0 if motif.upper() in aso_sequence.upper() else 0.0


def add_rnaseh1_motif_features(
        df: pd.DataFrame,
        out_prefix: str = "RNaseH1_"
) -> tuple[pd.DataFrame, list[str]]:
    """
    Extracts the longest DNA gap and checks for specific motifs known to
    enhance (Py-rich) or block (G-quadruplex) RNase H1 efficiency.

    Excludes pure toxicity markers (like TGC).
    """
    # Define motifs: Potency (High Efficacy) and Structural Blocks (Inefficiency)
    rnase_h_motifs = {
        # POSITIVE (High Potency) - Validated in Kiełpiński et al.
        'Potency_TCCC': 'TCCC',
        'Potency_TTCC': 'TTCC',
        'Potency_TCTC': 'TCTC',

        # NEGATIVE - inefficient
        'Inefficacy_GGGG': 'GGGG',
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
        df[col_name] = gap_sequences.apply(lambda x: check_motif_presence(x, motif_seq))

    return df, feature_cols