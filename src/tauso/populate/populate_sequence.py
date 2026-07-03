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
    ("seq_internal_fold", internal_fold),
    ("seq_internal_rna_fold", internal_rna_fold),
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


TERMINAL_5P = 6  # one-hot the first 6 bases from the 5' end: ohe_pos0.._pos5
TERMINAL_3P = 6  # and the first 6 bases from the 3' end:     ohe_3p0.._3p5


def populate_sequence_one_hot_encoded(
    df: pd.DataFrame, n5: int = TERMINAL_5P, n3: int = TERMINAL_3P, cpus: int = 1
) -> Tuple[pd.DataFrame, list[str]]:
    """One-hot encodes the ASO's TERMINAL bases only.

    The first `n5` bases from the 5' end become ohe_pos{0..n5-1}_{A,C,G,T}, and the first `n3` bases from
    the 3' end become ohe_3p{0..n3-1}_{A,C,G,T} (ohe_3p0 = the 3'-terminal base). U is encoded on the T
    channel; an unknown base or an out-of-range position (ASO shorter than the block) is all-zero. Anchoring
    the second block at the 3' terminus keeps it aligned across ASO lengths (16-mer cEt vs 20-mer MOE
    gapmers), which a single 5'-anchored absolute grid cannot do.

    Returns the updated DataFrame and the list of new feature names.
    """
    start_time = time.time()
    seqs = df[ASO_SEQUENCE].fillna("").astype(str).str.upper()

    def onehot(chars, base):
        hit = chars == base
        if base == "T":  # RNA U shares the T channel
            hit = hit | (chars == "U")
        return hit.to_numpy().astype("int8")

    columns, feature_names = {}, []
    for i in range(n5):  # 5' end: absolute positions 0..n5-1
        chars = seqs.str.slice(i, i + 1)  # empty when the ASO is shorter than i+1
        for nuc in ("A", "C", "G", "T"):
            name = f"ohe_pos{i}_{nuc}"
            columns[name] = onehot(chars, nuc)
            feature_names.append(name)
    for o in range(n3):  # 3' end: o-th base in from the 3' terminus
        chars = seqs.str[-1 - o]  # NaN when the ASO is shorter than o+1
        for nuc in ("A", "C", "G", "T"):
            name = f"ohe_3p{o}_{nuc}"
            columns[name] = onehot(chars, nuc)
            feature_names.append(name)

    encoded_df = pd.DataFrame(columns, index=df.index)
    existing_cols = [c for c in feature_names if c in df.columns]
    if existing_cols:
        df = df.drop(columns=existing_cols)
    df = pd.concat([df, encoded_df], axis=1)

    logger.debug("Finished terminal one-hot | Added %d features | Time: %.2fs",
                 len(feature_names), time.time() - start_time)
    return df, feature_names
