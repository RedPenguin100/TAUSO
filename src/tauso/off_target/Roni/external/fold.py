import ViennaRNA as RNA
import numpy as np

import os
import re
import random
import subprocess
import math

from Bio.Seq import Seq
from typing import Dict, List
from enum import Enum

from src.tauso.off_target.Roni.external.my_consts import RISEARCH1_BINARY_PATH, TMP_PATH ,RISEARCH1_PATH


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


class Interaction(Enum):
    RNA_RNA = "RNA_RNA"
    DNA_RNA_NO_WOBBLE = "DNA_RNA_NO_WOBBLE"
    DNA_DNA = "DNA_DNA"
    MODIFIED = "MODIFIED"


def dump_target_file(target_filename: str, name_to_sequence: Dict[str, str]):
    tmp_path = TMP_PATH / target_filename
    TMP_PATH.mkdir(exist_ok=True)
    with open(tmp_path, "w") as f:
        for name, sequence in name_to_sequence.items():
            f.write(f">{name}\n{sequence}\n")
    return tmp_path


def get_trigger_mfe_scores_by_risearch(trigger: str, name_to_sequence: Dict[str, str],
                                       interaction_type: Interaction = Interaction.DNA_RNA_NO_WOBBLE,
                                       minimum_score: int = 900, neighborhood: int = 0, parsing_type=None,
                                       target_file_cache=None,e_modified = (249,100)) -> str:
    if not RISEARCH1_BINARY_PATH.is_file():
        raise FileNotFoundError(
            f'RIsearch binary is not found at {RISEARCH1_BINARY_PATH}. Please read the "RIsearch" part of the project readme.')
    # used to dump cached files
    TMP_PATH.mkdir(exist_ok=True)

    hash = random.getrandbits(64)

    if target_file_cache is None:
        target_filename = f"target-{hash}.fa"
        target_path = dump_target_file(target_filename, name_to_sequence)
    else:
        target_filename = target_file_cache
        target_path = target_filename

    query_filename = f"query-{hash}.fa"
    query_path = TMP_PATH / query_filename

    with open(query_path, "w") as f:
        f.write(f">trigger\n{Seq(trigger).reverse_complement_rna()}\n")

    if interaction_type == Interaction.DNA_RNA_NO_WOBBLE:
        m = 'su95_noGU'
    elif interaction_type == Interaction.RNA_RNA:
        m = 't04'
    elif interaction_type == Interaction.MODIFIED:
        e1 , e2  = e_modified
        m = f'modified {e1} {e2}'  
    else:
        raise ValueError(f"Unsupported interaction type: {interaction_type}")

    try:
        args = [RISEARCH1_BINARY_PATH, "-q", str(query_path), "-t", str(target_path), "-s", f"{minimum_score}",
                "-d", "30", "-m", m, '-n', f"{neighborhood}"]
        if parsing_type is not None:
            args.append(f'-p{parsing_type}')

        result = subprocess.check_output(args, universal_newlines=True, text=True, cwd=str(RISEARCH1_PATH) 
                                         )
    finally:
        if target_file_cache is None:
            os.remove(target_path)
        query_path.unlink()
    # logger.info(f"finished evaluating a single trigger: {trigger}")
    return result
