import logging
import time
from typing import Iterable, Optional, Tuple

from Bio.SeqUtils import gc_fraction

logger = logging.getLogger(__name__)

from ..data.consts import ASO_SEQUENCE
from ..features.sequence.seq_features import *
from ..parallel_utils import make_apply_fn
from ..timer import Timer

logger = logging.getLogger(__name__)


FEATURE_SPECS: list[tuple[str, callable]] = [
    # Terminal Clamps
    # Returns 1 if 5' end is G/C, else 0
    ("seq_5_prime_clamp", lambda x: 1.0 if x and x[0].upper() in "GC" else 0.0),
    # Returns 1 if 3' end is G/C, else 0
    ("seq_3_prime_clamp", lambda x: 1.0 if x and x[-1].upper() in "GC" else 0.0),
    # ASO sequence energy
    ("seq_self_energy", self_energy),
    ("seq_internal_fold_dna", internal_fold),
    ("seq_internal_fold_rna", internal_fold_rna),
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


def calc_feature(df: pd.DataFrame, col_name: str, func, seq_series, cpus: int = 1, verbose=False) -> None:
    """
    Computes a feature with timing logs. Uses parallel_apply if available.

    ``seq_series`` is the ASO sequence normalized to uppercase DNA, so every feature is
    robust to lowercase or RNA (U) input.
    """
    start_time = time.time()
    logger.debug("Starting %s...", col_name)
    df[col_name] = make_apply_fn(seq_series, n_jobs=cpus)(func)
    duration = time.time() - start_time
    logger.debug("Finished %s | Time: %.2fs", col_name, duration)


def populate_sequence_features(
    df,
    features: Optional[Iterable[str]] = None,
    cpus: int = 1,
) -> Tuple:
    available = {name: fn for name, fn in FEATURE_SPECS}
    feature_names = list(features) if features is not None else [name for name, _ in FEATURE_SPECS]

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

    for name in feature_names:
        # Wrap each feature calculation in the Timer context manager
        with Timer(name):
            calc_feature(df, name, available[name], seq_series, cpus=cpus)

    return df, feature_names


import time
from typing import Optional, Tuple

import pandas as pd


def populate_sequence_one_hot_encoded(
    df: pd.DataFrame, max_len: Optional[int] = None, cpus: int = 1
) -> Tuple[pd.DataFrame, list[str]]:
    """
    One-hot encodes ASO sequences with zero-padding to handle varying lengths.
    Flattens the output into tabular columns (e.g., ohe_pos0_A, ohe_pos0_C...).

    Returns:
        Tuple containing the updated DataFrame and a list of the new feature names.
    """
    start_time = time.time()

    logger.debug("Starting One-Hot Encoding with Zero-Padding...")

    # 1. Determine max length for padding (if not manually provided)
    if max_len is None:
        max_len = int(df[ASO_SEQUENCE].str.len().max())
        logger.debug("Determined max_len automatically: %d", max_len)

    # 2. Map nucleotides to lists (includes U for RNA compatibility)
    # Unknowns or Ns will map to [0, 0, 0, 0]
    nuc_map = {
        "A": [1, 0, 0, 0],
        "C": [0, 1, 0, 0],
        "G": [0, 0, 1, 0],
        "T": [0, 0, 0, 1],
        "U": [0, 0, 0, 1],
    }
    pad_vec = [0, 0, 0, 0]

    def encode_and_pad(seq: str) -> list[int]:
        if not isinstance(seq, str):
            seq = ""

        encoded = []
        # Encode up to max_len (truncates if sequence exceeds max_len)
        for nuc in seq[:max_len]:
            encoded.extend(nuc_map.get(nuc.upper(), pad_vec))

        # Post-pad with zeros if sequence is shorter than max_len
        pad_length = max_len - len(seq)
        if pad_length > 0:
            encoded.extend(pad_vec * pad_length)

        return encoded

    # 3. Apply function
    encoded_series = make_apply_fn(df[ASO_SEQUENCE], n_jobs=cpus)(encode_and_pad)

    # 4. Generate the new feature names
    feature_names = []
    for pos in range(max_len):
        for nuc in ["A", "C", "G", "T"]:
            feature_names.append(f"ohe_pos{pos}_{nuc}")

    # 5. Expand the lists into DataFrame columns
    # (Doing this via pd.DataFrame conversion is significantly faster than looping column-wise)
    encoded_df = pd.DataFrame(encoded_series.tolist(), columns=feature_names, index=df.index)

    # 6. Safety check: Drop existing OHE columns if re-running to avoid duplication
    existing_cols = [c for c in feature_names if c in df.columns]
    if existing_cols:
        df = df.drop(columns=existing_cols)

    # Concat the new features back to the main dataframe
    df = pd.concat([df, encoded_df], axis=1)

    duration = time.time() - start_time
    logger.debug("Finished One-Hot Encoding | Added %d features | Time: %.2fs", len(feature_names), duration)

    return df, feature_names
