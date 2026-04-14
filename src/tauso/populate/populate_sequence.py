import time
from typing import Iterable, Optional, Tuple

from Bio.SeqUtils import gc_fraction

from ..data.consts import SEQUENCE
from ..features.sequence.seq_features import *
from ..timer import Timer

FEATURE_SPECS: list[tuple[str, callable]] = [
    # Terminal Clamps
    # Returns 1 if 5' end is G/C, else 0
    ("Sequence_5_prime_clamp", lambda x: 1.0 if x and x[0].upper() in "GC" else 0.0),
    # Returns 1 if 3' end is G/C, else 0
    ("Sequence_3_prime_clamp", lambda x: 1.0 if x and x[-1].upper() in "GC" else 0.0),
    # ASO sequence energy
    ("Sequence_self_energy", self_energy),
    ("Sequence_internal_fold", internal_fold),
    # Basic Composition Features
    ("Sequence_purine_content", purine_content),
    ("Sequence_gc_content", gc_fraction),
    ("Sequence_ggg_counts", count_g_runs),
    ("Sequence_a_count", lambda x: x.count("A")),
    ("Sequence_c_count", lambda x: x.count("C")),
    ("Sequence_g_count", lambda x: x.count("G")),
    ("Sequence_t_count", lambda x: x.count("T") + x.count("U")),  # In case the sequence is with U for some reason.
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
    ("Sequence_ta_dinucleotide_fraction", ta_dinucleotide_fraction),
    ("Sequence_poly_pyrimidine_stretch", poly_pyrimidine_stretch),
    ("Sequence_gc_block_length", gc_block_length),
    ("Sequence_at_rich_region_score", at_rich_region_score),
    ("Sequence_gc_content_3prime_end", gc_content_3prime_end),
    ("Sequence_homooligo_count", homooligo_count),
    ("Sequence_absolute_terminal_gc", absolute_terminal_gc),
    ("Sequence_keto_amino_skew", keto_amino_skew),
    ("Sequence_ry_transition_fraction", ry_transition_fraction),
]


def calc_feature(df: pd.DataFrame, col_name: str, func, cpus: int = 1, verbose=False) -> None:
    """
    Computes a feature with timing logs. Uses parallel_apply if available.
    """
    start_time = time.time()

    if verbose:
        print(f"► Starting {col_name}...")

    series = df[SEQUENCE]
    try:
        from pandarallel import pandarallel

        pandarallel.initialize(progress_bar=verbose, verbose=0, nb_workers=cpus)
        df[col_name] = series.parallel_apply(func)
    except Exception:
        df[col_name] = series.apply(func)

    if verbose:
        duration = time.time() - start_time
        print(f"✔ Finished {col_name} | Time: {duration:.2f}s\n")


def populate_sequence_features(
    df,
    features: Optional[Iterable[str]] = None,
    cpus: int = 1,
) -> Tuple:
    available = {name: fn for name, fn in FEATURE_SPECS}
    feature_names = list(features) if features is not None else [name for name, _ in FEATURE_SPECS]

    for name in feature_names:
        # Wrap each feature calculation in the Timer context manager
        with Timer(name):
            calc_feature(df, name, available[name], cpus=cpus)

    return df, feature_names


def populate_sequence_one_hot_encoded(
    df: pd.DataFrame,
    max_len: Optional[int] = None,
    cpus: int = 1,
    verbose: bool = False,
) -> Tuple[pd.DataFrame, list[str]]:
    """
    One-hot encodes ASO sequences with zero-padding to handle varying lengths.
    Flattens the output into tabular columns (e.g., OHE_pos0_A, OHE_pos0_C...).

    Returns:
        Tuple containing the updated DataFrame and a list of the new feature names.
    """
    start_time = time.time()

    if verbose:
        print(f"► Starting One-Hot Encoding with Zero-Padding...")

    # 1. Determine max length for padding (if not manually provided)
    if max_len is None:
        max_len = int(df[SEQUENCE].str.len().max())
        if verbose:
            print(f"  ↳ Determined max_len automatically: {max_len}")

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

    # 3. Apply function (mirroring your pandarallel structure)
    try:
        from pandarallel import pandarallel

        # Only initialize if not already done, though harmless if repeated
        pandarallel.initialize(progress_bar=verbose, verbose=0, nb_workers=cpus)
        encoded_series = df[SEQUENCE].parallel_apply(encode_and_pad)
    except Exception:
        encoded_series = df[SEQUENCE].apply(encode_and_pad)

    # 4. Generate the new feature names
    feature_names = []
    for pos in range(max_len):
        for nuc in ["A", "C", "G", "T"]:
            feature_names.append(f"OHE_pos{pos}_{nuc}")

    # 5. Expand the lists into DataFrame columns
    # (Doing this via pd.DataFrame conversion is significantly faster than looping column-wise)
    encoded_df = pd.DataFrame(encoded_series.tolist(), columns=feature_names, index=df.index)

    # 6. Safety check: Drop existing OHE columns if re-running to avoid duplication
    existing_cols = [c for c in feature_names if c in df.columns]
    if existing_cols:
        df = df.drop(columns=existing_cols)

    # Concat the new features back to the main dataframe
    df = pd.concat([df, encoded_df], axis=1)

    if verbose:
        duration = time.time() - start_time
        print(f"✔ Finished One-Hot Encoding | Added {len(feature_names)} features | Time: {duration:.2f}s\n")

    return df, feature_names
