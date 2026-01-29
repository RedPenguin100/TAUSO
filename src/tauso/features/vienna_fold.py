import ViennaRNA as RNA
import numpy as np

import math

def calculate_avg_mfe_over_sense_region(sequence, sense_start, sense_length, flank_size=120, window_size=120, step=1):
    sequence = str(sequence).upper().replace('T', 'U')
    sequence_length = len(sequence)
    energy_values = np.zeros(sequence_length)
    counts = np.zeros(sequence_length)

    for i in range(0, sequence_length - window_size + 1, step):
        subseq = sequence[i:i + window_size]
        fc = RNA.fold_compound(subseq)
        _, mfe = fc.mfe()
        mfe_per_nt = mfe / window_size

        for j in range(i, i + window_size):
            energy_values[j] += mfe_per_nt
            counts[j] += 1

    counts[counts == 0] = 1
    avg_energies = energy_values / counts

    flank_start = max(0, sense_start - flank_size)
    sense_start_in_flank = sense_start - flank_start
    sense_end_in_flank = sense_start_in_flank + sense_length

    if 0 <= sense_start_in_flank < sequence_length and sense_end_in_flank <= sequence_length:
        return np.mean(avg_energies[sense_start_in_flank:sense_end_in_flank])
    else:
        return np.nan



def get_weighted_energy(target_start, l, step_size, energies, window_size):
    """
    Calculate average energy for a target region by finding which sliding windows
    overlap each position and averaging their energies.
    """
    if l <= 0 or len(energies) == 0:
        return 0.0

    num_windows = len(energies)
    position_energies = np.zeros(l, dtype=np.float64)

    for position in range(target_start, target_start + l):
        # Find all windows that overlap this position
        # A window at index k covers positions from k*step_size to k*step_size + window_size - 1

        # First window that could overlap: its end must reach this position
        first_window = max(0, math.ceil((position - window_size + 1) / step_size))

        # Last window that could overlap: its start must be at or before this position
        last_window = min(num_windows - 1, position // step_size)

        if first_window > last_window:
            # No windows fully overlap - use the closest window
            closest_window = np.clip(position // step_size, 0, num_windows - 1)
            energy_value = float(energies[closest_window])
        else:
            # Average all overlapping windows
            energy_value = float(np.mean(energies[first_window:last_window + 1], dtype=np.float64))

        position_energies[position - target_start] = energy_value

    return float(np.mean(position_energies, dtype=np.float64))


def calculate_energies(target_seq, step_size, window_size):
    L = len(target_seq)
    if L < window_size: return np.empty(0, dtype=np.float64)

    starts = list(range(0, L - window_size + 1, step_size))
    last_needed = L - window_size
    if starts[-1] != last_needed:
        starts.append(last_needed)  # add [L-window_size, L-1] window

    energies = np.empty(len(starts), dtype=np.float64)
    for k, i in enumerate(starts):
        _, mfe = RNA.fold(target_seq[i:i + window_size])
        energies[k] = float(mfe)
    return energies
