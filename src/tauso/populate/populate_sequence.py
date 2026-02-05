from typing import Iterable, Optional, Tuple

from Bio.SeqUtils import gc_fraction

from ..features.sequence.seq_features import *


FEATURE_SPECS: list[tuple[str, callable]] = [
    ("Sequence_self_energy", self_energy),
    ("Sequence_internal_fold", internal_fold),

    # Basic Composition Features
    ("Sequence_purine_content", purine_content),
    ("Sequence_gc_content", gc_fraction),
    ("Sequence_ggg_counts", count_g_runs),

    # Palindromes and Entropy
    ("Sequence_4_palindromic", lambda x: palindromic_fraction(x, 4)),
    ("Sequence_6_palindromic", lambda x: palindromic_fraction(x, 6)),
    ("Sequence_entropy", seq_entropy),
    ("Sequence_dinucleotide_entropy", dinucleotide_entropy),
    ("Sequence_nucleotide_diversity", nucleotide_diversity),

    # Structure and Energy (Hairpins/Motifs)
    ("Sequence_hairpin_score", hairpin_score),
    ("Sequence_hairpin_dG_energy", hairpin_dG_energy),
    ("Sequence_hairpin_tm", hairpin_tm),
    ("Sequence_toxic_motif_count", toxic_motif_count),

    # Repeats and Skews
    ("Sequence_tandem_repeats_score", tandem_repeats_score),
    ("Sequence_dispersed_repeats_score", dispersed_repeats_score),
    ("Sequence_gc_skew", gc_skew),
    ("Sequence_gc_skew_ends", gc_skew_ends),
    ("Sequence_at_skew", at_skew),

    # Advanced Region Scores
    ("Sequence_flexible_dinucleotide_fraction", flexible_dinucleotide_fraction),
    ("Sequence_stop_codon_count", stop_codon_count),
    ("Sequence_cg_dinucleotide_fraction", cg_dinucleotide_fraction),
    ("Sequence_poly_pyrimidine_stretch", poly_pyrimidine_stretch),
    ("Sequence_gc_block_length", gc_block_length),
    ("Sequence_at_rich_region_score", at_rich_region_score),
    ("Sequence_gc_content_3prime_end", gc_content_3prime_end),
    ("Sequence_homooligo_count", homooligo_count),
]


def populate_sequence_features(
    df,
    features: Optional[Iterable[str]] = None,
    input_col: str = "Sequence",
    cpus: int = 1,
) -> Tuple:
    """
    Populate sequence feature columns on the provided DataFrame.

    Args:
        df: Input DataFrame to populate in-place.
        features: Optional iterable of feature names to compute. Defaults to all features.
        input_col: Column name containing sequences.
        cpus: Number of CPUs to use for parallel feature computation (default=1).

    Returns:
        (df, feature_names) tuple.
    """
    available = {name: fn for name, fn in FEATURE_SPECS}
    if features is None:
        feature_names = [name for name, _ in FEATURE_SPECS]
    else:
        feature_names = list(features)
        unknown = [name for name in feature_names if name not in available]
        if unknown:
            raise ValueError(f"Unknown sequence features: {unknown}")

    for feature_name in feature_names:
        calc_feature(df, feature_name, available[feature_name], input_col=input_col, cpus=cpus)

    return df, feature_names
