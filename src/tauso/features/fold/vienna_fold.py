from functools import lru_cache

import numpy as np
import ViennaRNA as RNA


@lru_cache(maxsize=10000)
def get_cached_mfe(seq):
    return RNA.fold(seq)[1]


def _per_position_avg_energies(sequence, window_size, start, stop, step):
    """One sliding-window sweep. Returns per-position avg mfe across `sequence`.

    Same inner loop as the single-step reference: fold each window in
    range(start, stop, step), spread mfe/window_size over its positions, then
    average per position. Positions not covered by any window are NaN.
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

    Equivalent to calling the single-step reference once per `step` in `steps`, but faster:
      1. We only sweep windows whose `[i, i+window_size)` overlaps the sense region.
         Windows outside contribute to positions never read by the final
         `nanmean(avg[sense_start:sense_end])`, so skipping them is identity-preserving.
      2. Different steps share fold calls at overlapping positions via the module-level
         `get_cached_mfe` lru_cache (within a single call the unique-window count is
         well under cache size, so this is effectively free).
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
