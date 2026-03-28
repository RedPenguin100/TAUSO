from functools import lru_cache

import numpy as np
import ViennaRNA as RNA


# --- NEW: Helper with Caching (Strategy 4) ---
# This makes repetitive calls (like in your test) instant.
# We use RNA.fold instead of fold_compound to save memory/init time (Strategy 2).
@lru_cache(maxsize=10000)
def get_cached_mfe(seq):
    return RNA.fold(seq)[1]


def calculate_avg_mfe_over_sense_region(sequence, sense_start_in_flank, sense_length, window_size=120, step=1):
    sequence = str(sequence).upper().replace("T", "U")
    sequence_length = len(sequence)
    energy_values = np.zeros(sequence_length)
    counts = np.zeros(sequence_length)

    # --- CHANGED: Standard Loop using Cache ---
    # We reverted to the loop (to fix the crash) but use the cached helper.
    for i in range(0, sequence_length - window_size + 1, step):
        subseq = sequence[i: i + window_size]

        # This call is now cached and overhead-free
        mfe = get_cached_mfe(subseq)

        # You were ALREADY doing Strategy 3 (Vectorized Slicing) here. Good job!
        mfe_per_nt = mfe / window_size
        energy_values[i: i + window_size] += mfe_per_nt
        counts[i: i + window_size] += 1
    # ------------------------------------------

    avg_energies = np.divide(
        energy_values,
        counts,
        out=np.full_like(energy_values, np.nan),
        where=counts != 0,
    )

    sense_end_in_flank = sense_start_in_flank + sense_length
    if (
            0 <= sense_start_in_flank < sequence_length
            and sense_end_in_flank <= sequence_length
    ):
        return np.nanmean(avg_energies[sense_start_in_flank:sense_end_in_flank])
    else:
        return np.nan
