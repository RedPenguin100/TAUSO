import ViennaRNA as RNA
import numpy as np

import re
import math

from typing import List

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


'''
note: I used hashing here to be able to run multiple sequences (triggers in my case) in parallel (for multiple triggers) without conflicts with the files,
you can remove it if it's not an issue for you 
note: in my case I first preformed reverse_complement_rna to the trigger because my goal was to find sequences that might undesirably bind to the
trigger binding site of the toehold, adjust it for your particular case
note: play with the d and s parameters I passed here, and with other parameters that might be relevant for your usecase
'''


def _parse_mfe_scores_2(result):
    if not result:
        return [[]]

    lines = result.split('\n')
    target_to_energies = dict()
    for line in lines:
        if line == '':
            continue
        line_parts = line.split('\t')
        target_name = line_parts[3]
        target_energy = line_parts[-1]
        if line_parts[3] in target_to_energies:
            target_to_energies[target_name].append(float(target_energy))
        else:
            target_to_energies[target_name] = [float(target_energy)]

    # ignore the keys at this point
    return list(target_to_energies.values())


# in my case I was only interested in the energy scores and didn't care about the actual sequences, this is the parsing I used in case it helps you
def get_mfe_scores(result: str, parsing_type=None) -> List[List[float]]:
    mfe_results = []

    if parsing_type is None:
        for gene_result in result.split("\n\nquery trigger")[1:]:
            stripped_result = gene_result.strip()
            regex_results = re.findall("Free energy \[kcal/mol\]: [0-9-.]+ ", stripped_result)
            mfe_results.append(
                [float(regex_result.replace('Free energy [kcal/mol]: ', '').strip()) for regex_result in regex_results])
    elif parsing_type == '2':
        return _parse_mfe_scores_2(result)

    return mfe_results

