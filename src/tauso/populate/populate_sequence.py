import logging
import time
from typing import Iterable, Optional, Tuple

import pandas as pd
from Bio.SeqUtils import gc_fraction

logger = logging.getLogger(__name__)

from ..data.consts import ASO_SEQUENCE
from ..features.sequence.seq_features import *
from ..parallel_utils import make_apply_fn
from .feature_runner import FeatureSpec, compute_features

FEATURE_SPECS: list[FeatureSpec] = [
    # Terminal Clamps
    # Returns 1 if 5' end is G/C, else 0
    ("seq_5_prime_clamp", lambda x: 1.0 if x and x[0].upper() in "GC" else 0.0),
    # Returns 1 if 3' end is G/C, else 0
    ("seq_3_prime_clamp", lambda x: 1.0 if x and x[-1].upper() in "GC" else 0.0),
    # ASO sequence energy
    ("seq_self_energy", self_energy),
    ("seq_internal_fold_rna", internal_fold_rna),
    ("sense_internal_fold", sense_internal_fold),
    # Basic Composition Features
    ("seq_purine_content", purine_content),
    ("seq_gc_content", gc_fraction),
    ("seq_a_count", lambda x: x.upper().count("A")),
    ("seq_c_count", lambda x: x.upper().count("C")),
    ("seq_g_count", lambda x: x.upper().count("G")),
    (
        "seq_t_count",
        lambda x: x.upper().count("T") + x.upper().count("U"),
    ),  # In case the sequence is with U for some reason.
    # Palindromes and Entropy
    ("seq_4_palindromic", lambda x: palindromic_fraction(x, 4)),
    ("seq_6_palindromic", lambda x: palindromic_fraction(x, 6)),
    ("seq_entropy", seq_entropy),
    ("seq_dinucleotide_entropy", dinucleotide_entropy),
    ("seq_dinucleotide_diversity", dinucleotide_diversity),
    # Structure and Energy (Hairpins/Motifs)
    ("seq_hairpin_score", hairpin_score),
    ("seq_hairpin_dG_energy", hairpin_dG_energy),
    ("seq_hairpin_tm", hairpin_tm),
    # Repeats and Skews
    ("seq_tandem_repeats_score", tandem_repeats_score),
    ("seq_dispersed_repeats_score", dispersed_repeats_score),
    ("seq_gc_skew", gc_skew),
    ("seq_gc_skew_ends", gc_skew_ends),
    ("seq_at_skew", at_skew),
    # Advanced Region Scores
    ("seq_flexible_dinucleotide_fraction", flexible_dinucleotide_fraction),
    ("seq_cg_dinucleotide_fraction", cg_dinucleotide_fraction),
    ("seq_ta_dinucleotide_fraction", ta_dinucleotide_fraction),
    ("seq_poly_pyrimidine_run_count", poly_pyrimidine_run_count),
    ("seq_poly_pyrimidine_longest", poly_pyrimidine_longest),
    ("seq_gc_longest", gc_longest),
    ("seq_at_rich_run_count", at_rich_run_count),
    ("seq_at_rich_longest", at_rich_longest),
    ("seq_gc_content_3prime_end", gc_content_3prime_end),
    ("seq_homooligo_count", homooligo_count),
    ("seq_terminal_gc_fraction", terminal_gc_fraction),
    ("seq_keto_amino_skew", keto_amino_skew),
    ("seq_ry_transition_fraction", ry_transition_fraction),
]


def populate_sequence_features(
    df,
    features: Optional[Iterable[str]] = None,
    cpus: int = 1,
) -> Tuple:
    # Normalize to uppercase DNA once so every feature is robust to lowercase or RNA (U) input.
    seq_series = df[ASO_SEQUENCE].str.upper().str.replace("U", "T", regex=False)

    # Disallow anything but A/C/G/T (e.g. N, gaps, IUPAC ambiguity codes). A feature computed on an
    # undefined base is silently wrong, which is unacceptable for a design tool — fail loudly instead.
    invalid = seq_series[seq_series.str.contains(r"[^ACGT]", regex=True, na=True)]
    if len(invalid) > 0:
        raise ValueError(
            f"{len(invalid)} ASO sequence(s) contain non-ACGT characters after upper/U->T normalization "
            f"(only A/C/G/T/U accepted). Examples: {invalid.head(5).tolist()}"
        )

    return compute_features(
        df,
        FEATURE_SPECS,
        seq_series,
        lambda apply_fn, func: apply_fn(func),
        features=features,
        cpus=cpus,
    )


# Nucleotides one-hot-encoded inward from each terminus of the ASO. Both ends are encoded in
# their own frame -- 5' positions from the start, 3' positions inward from the last base -- so
# the terminal signal (5' seed, 3' RNase-H gap edge) stays aligned across variable-length
# gapmers. Each position is assigned to whichever half of the sequence it falls in; positions
# past a sequence's midpoint (or its end) are left NaN, so the two windows never double-encode
# the middle of a short ASO.
TERMINAL_OHE_N = 10


def one_hot_feature_names(terminal_n: int = TERMINAL_OHE_N) -> list[str]:
    """Column names for terminal one-hot encoding: the first ``terminal_n`` nucleotides
    from the 5' end (``ohe_pos{i}_{nuc}``) followed by the first ``terminal_n`` from the
    3' end (``ohe_3p{i}_{nuc}``, counting inward from each ASO's own 3' terminus)."""
    names = []
    for pos in range(terminal_n):
        for nuc in ("A", "C", "G", "T"):
            names.append(f"ohe_pos{pos}_{nuc}")
    for pos in range(terminal_n):
        for nuc in ("A", "C", "G", "T"):
            names.append(f"ohe_3p{pos}_{nuc}")
    return names


def populate_sequence_one_hot_encoded(
    df: pd.DataFrame, terminal_n: int = TERMINAL_OHE_N, cpus: int = 1
) -> Tuple[pd.DataFrame, list[str]]:
    """
    One-hot encodes the terminal nucleotides of each ASO sequence.

    Encodes up to ``terminal_n`` nucleotides from the 5' end into ``ohe_pos{i}_{nuc}`` columns
    and up to ``terminal_n`` from the 3' end into ``ohe_3p{i}_{nuc}`` columns. The 3' window is
    length-aware: ``ohe_3p{i}`` is the nucleotide ``i`` positions inward from each sequence's
    own 3' terminus, so it aligns across variable-length gapmers. Each nucleotide is assigned to
    the terminus on its own half of the sequence; positions past the midpoint (or the end) of a
    short ASO are NaN, so the 5' and 3' windows never encode the same base twice.

    Returns:
        Tuple containing the updated DataFrame and a list of the new feature names.
    """
    start_time = time.time()

    logger.debug("Starting terminal One-Hot Encoding (%d nt per end)...", terminal_n)

    # Map nucleotides to lists (includes U for RNA compatibility).
    # Unknowns or Ns, and positions outside a terminus's own half of the sequence, map to NaN.
    nuc_map = {
        "A": [1, 0, 0, 0],
        "C": [0, 1, 0, 0],
        "G": [0, 0, 1, 0],
        "T": [0, 0, 0, 1],
        "U": [0, 0, 0, 1],
    }
    nan_vec = [float("nan")] * 4

    def encode_terminal(seq: str) -> list[float]:
        if not isinstance(seq, str):
            seq = ""
        seq = seq.upper()
        n = len(seq)

        encoded = []
        # 5' terminus: nucleotide i from the start, kept only while it lies in the 5' half.
        for i in range(terminal_n):
            encoded.extend(nuc_map.get(seq[i], nan_vec) if (i < n and 2 * i <= n - 1) else nan_vec)
        # 3' terminus: nucleotide i inward from the last base, kept only in the 3' half.
        for i in range(terminal_n):
            j = n - 1 - i
            encoded.extend(nuc_map.get(seq[j], nan_vec) if (j >= 0 and 2 * j > n - 1) else nan_vec)

        return encoded

    # Apply function
    encoded_series = make_apply_fn(df[ASO_SEQUENCE], n_jobs=cpus)(encode_terminal)

    feature_names = one_hot_feature_names(terminal_n)

    # Expand the lists into DataFrame columns
    # (Doing this via pd.DataFrame conversion is significantly faster than looping column-wise)
    encoded_df = pd.DataFrame(encoded_series.tolist(), columns=feature_names, index=df.index)

    # Safety check: drop existing OHE columns if re-running to avoid duplication
    existing_cols = [c for c in feature_names if c in df.columns]
    if existing_cols:
        df = df.drop(columns=existing_cols)

    # Concat the new features back to the main dataframe
    df = pd.concat([df, encoded_df], axis=1)

    duration = time.time() - start_time
    logger.debug("Finished terminal One-Hot Encoding | Added %d features | Time: %.2fs", len(feature_names), duration)

    return df, feature_names
