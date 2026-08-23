import numpy as np

from ...util import dna_to_rna
from ._vienna_internal import window_mfe


def _get_cached_mfe(sequence, start=0, window_size=None):
    """MFE of ``sequence[start:start + window_size]``, sharing one fold of `sequence`.

    The whole `sequence` is folded once and cached, and each window is read back out of
    those shared matrices rather than folded on its own.
    """
    if window_size is None:
        window_size = len(sequence) - start
    return window_mfe(sequence, start, window_size)


def _per_position_avg_energies(sequence, window_size, start, stop, step):
    """One sliding-window sweep over `sequence`; returns the per-position average MFE.

    Fold each window in range(start, stop, step), spread its mfe/window_size
    evenly over the window's positions, then average per position. Positions no
    window covers are NaN.
    """
    seq_len = len(sequence)
    energy_values = np.zeros(seq_len)
    counts = np.zeros(seq_len)
    for i in range(start, stop, step):
        # get_catched_mfe is equivalent to RNA.fold(sequence[i:i+window_size])[1], but faster.
        # only fold the entire sequence once, and then fetch folded items.
        mfe_per_nt = _get_cached_mfe(sequence, i, window_size) / window_size
        energy_values[i : i + window_size] += mfe_per_nt
        counts[i : i + window_size] += 1
    return np.divide(
        energy_values,
        counts,
        out=np.full_like(energy_values, np.nan),
        where=counts != 0,
    )


def calculate_avg_mfe_per_step(sequence, sense_start_in_flank, sense_length, window_size, steps):
    """Sliding-window MFE averaged over the sense region, for several `step` values at once.

    The plain version of each step sweeps windows from position 0 across the whole
    sequence and nanmeans the sense-region positions. Two shortcuts give the same
    numbers more cheaply:
      1. Only sweep windows overlapping the sense region — outside windows touch
         only positions the final mean never reads.
      2. Overlapping windows share a single fold of `sequence`, each reading its own
         energy back out of the shared matrices.
    """
    sequence = dna_to_rna(sequence)
    seq_len = len(sequence)
    sense_end = sense_start_in_flank + sense_length

    if not (0 <= sense_start_in_flank < seq_len and sense_end <= seq_len):
        return {step: np.nan for step in steps}

    # A window at position i overlaps the sense region iff
    # i+window-1 >= sense_start AND i <= sense_end-1.
    min_i = max(0, sense_start_in_flank - window_size + 1)
    max_i = min(seq_len - window_size, sense_end - 1)

    out = {}
    for step in steps:
        # First multiple of `step` >= min_i. Same grid as starting at 0, just filtered.
        start_i = ((min_i + step - 1) // step) * step
        avg_energies = _per_position_avg_energies(sequence, window_size, start_i, max_i + 1, step)
        out[step] = np.nanmean(avg_energies[sense_start_in_flank:sense_end])
    return out
