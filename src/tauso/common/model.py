from typing import Optional

import pandas as pd

from tauso.data.consts import *
from tauso.genome.read_human_genome import get_locus_to_data_dict
from tauso.util import get_antisense

MOD_TYPE_DICT = {'moe': 'MMMMMddddddddddMMMMM', 'lna': 'LLLddddddddddLLL'}

def get_target_sequence(
        target_name: str,
        custom_sequence: Optional[str] = None
) -> str:
    """
    Retrieves the mRNA sequence either from the custom input or
    by querying the local genome database.
    """
    if custom_sequence:
        return custom_sequence.upper().replace('T', 'U')

    # Use the optimized subset loader from read_human_genome.py
    print(f"Fetching sequence for {target_name}...")
    locus_data = get_locus_to_data_dict(include_introns=False, gene_subset=[target_name])

    if target_name not in locus_data:
        raise ValueError(f"Gene {target_name} not found in genome database.")

    # Extract the full mRNA sequence from the LocusInfo object
    return locus_data[target_name].full_mrna


DEFAULT_ASO_LENGTH = 20


def get_init_df(target_mrna: str, mod_type: str = 'moe') -> pd.DataFrame:
    """
    Generates ASO candidates by tiling the mRNA.

    Changes from legacy:
    - Removed 'end' parameter.
    - uses len(target_mrna) to calculate 'sense_start_from_end'.
    """
    # Determine ASO length (s) based on chemistry
    if mod_type in MOD_TYPE_DICT:
        s = len(MOD_TYPE_DICT[mod_type])
    else:
        s = 20  # Default fallback

    target_len = len(target_mrna)

    candidates = []
    sense_starts = []
    sense_lengths = []
    sense_starts_from_end = []

    # Deduplication set
    set_candidates = set()

    # Tile the sequence
    for i in range(0, target_len - (s - 1)):
        target = target_mrna[i: i + s]

        # Deduplicate based on the SENSE target sequence
        if target in set_candidates:
            continue
        set_candidates.add(target)

        # Store Data
        candidates.append(get_antisense(str(target)))
        sense_starts.append(i)
        sense_lengths.append(s)

        # distance from 3' end (assuming full_mrna is the context)
        # Note: Index 0 is 5' end. Distance from 3' end decreases as i increases.
        sense_starts_from_end.append(target_len - i)

    df = pd.DataFrame({
        SEQUENCE: candidates,
        SENSE_START: sense_starts,
        SENSE_LENGTH: sense_lengths,
        "sense_start_from_end": sense_starts_from_end
    })

    return df