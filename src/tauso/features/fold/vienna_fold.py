from functools import lru_cache

import numpy as np
import ViennaRNA as RNA


@lru_cache(maxsize=10000)
def get_cached_mfe(seq):
    return RNA.fold(seq)[1]


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
        mfe_per_nt = get_cached_mfe(sequence[i : i + window_size]) / window_size
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
      2. Overlapping windows reuse folds across steps via the module-level
         `get_cached_mfe` cache (unique windows per call stay well under its size).
    """
    sequence = str(sequence).upper().replace("T", "U")
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
