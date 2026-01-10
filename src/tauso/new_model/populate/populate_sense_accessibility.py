import numpy as np
import pandas as pd
from ..consts_dataframe import SENSE_START
from ...features.feature_names import SENSE_LENGTH
from ...features.rna_access.access_calculator import get_cache, get_sense_with_flanks
from ...features.rna_access.sense_accessibility import compute_sense_accessibility_value
from ...features.rna_access.access_calculator import AccessCalculator
SENSE_AVG_ACCESSIBILITY = 'sense_avg_accessibility'

# --- CONSTANTS ---
FLANK_SIZE = 120
ACCESS_SIZE = 13
SEED_SIZE = 13
SEED_SIZES = [SEED_SIZE * m for m in range(1, 4)]
ACCESS_WIN_SIZE = 80

def populate_sense_accessibility(aso_dataframe):

    if 'pre_mrna_sequence' not in aso_dataframe.columns:
        raise ValueError("DataFrame must contain 'pre_mrna_sequence' column used for extraction.")
    access_cache = get_cache(SEED_SIZES, access_size=ACCESS_SIZE)

    valid_rows = [
        (idx, row) for idx, row in aso_dataframe.iterrows()
        if row[SENSE_START] != -1
    ]

    for idx, row in valid_rows:
        sense_start = row[SENSE_START]
        sense_length = row[SENSE_LENGTH]
        full_mrna_seq = str(row['pre_mrna_sequence'])

        flanked_sense = get_sense_with_flanks(
            full_mrna_seq, sense_start, sense_length,
            flank_size=FLANK_SIZE
        )

        if not flanked_sense or len(flanked_sense) < ACCESS_SIZE:
            aso_dataframe.loc[idx, SENSE_AVG_ACCESSIBILITY] = 0
            continue

        avg_sense_access = compute_sense_accessibility_value(
            sense_start, sense_length, flank=flanked_sense, flank_size=FLANK_SIZE, access_win_size=ACCESS_WIN_SIZE,
            seed_sizes=SEED_SIZES, access_size=ACCESS_SIZE, cache=access_cache
        )
        aso_dataframe.loc[idx, SENSE_AVG_ACCESSIBILITY] = avg_sense_access

def calculate_sense_accessibility_batch(aso_dataframe, batch_size=1000):

    if 'pre_mrna_sequence' not in aso_dataframe.columns:
        raise ValueError("DataFrame must contain 'pre_mrna_sequence' column used for extraction.")
    access_cache = get_cache(SEED_SIZES, access_size=ACCESS_SIZE)

    valid_rows = [
        (idx, row) for idx, row in aso_dataframe.iterrows()
        if row[SENSE_START] != -1
    ]
    all_results = []
    # creating batches from the df
    for batch_start in range(0, len(valid_rows), batch_size):
        batch_rows = valid_rows[batch_start:batch_start + batch_size]

        rna_seqs = []
        sense_len_map = {}
        sense_start_map = {}

        for idx, row in batch_rows:
            full_mrna_seq = str(row['pre_mrna_sequence'])
            current_sense_start = row[SENSE_START]
            flank_start = max(0, current_sense_start - FLANK_SIZE)
            relative_start = current_sense_start - flank_start

            flanked_sense = get_sense_with_flanks(
                full_mrna_seq,
                current_sense_start,
                row[SENSE_LENGTH],
                flank_size=FLANK_SIZE
            )
            if not flanked_sense or len(flanked_sense) < ACCESS_SIZE:
                continue
            flanked_sense = flanked_sense.upper().replace('T', 'U')
            rna_seqs.append((str(idx), flanked_sense))
            sense_len_map[str(idx)] = row[SENSE_LENGTH]
            sense_start_map[str(idx)] = relative_start

        if not rna_seqs:
            continue

        # calculating the accessibility for a batch
        df = AccessCalculator.calc_new(
            rna_seqs,
            access_size=ACCESS_SIZE,
            min_gc=0,
            max_gc=100,
            gc_ranges=1,
            access_win_size=ACCESS_WIN_SIZE,
            access_seed_sizes=SEED_SIZES,
            cache=access_cache)


        df['pos_in_flank'] = df.groupby('rna_id').cumcount()
        df['sense_length'] = df['rna_id'].map(sense_len_map)
        df['rel_start'] = df['rna_id'].map(sense_start_map)
        mask = (
                (df['pos_in_flank'] >= df['rel_start']) &
                (df['pos_in_flank'] < df['rel_start'] + df['sense_length'])
        )
        batch_result = (
            df.loc[mask]
            .groupby('rna_id')['avg_access']
            .mean()
            .reset_index(name='access')
        )
        batch_result['rna_id'] = batch_result['rna_id'].astype(int)
        all_results.append(batch_result)

    if not all_results:
        return pd.DataFrame(columns=['rna_id', 'access'])

    return pd.concat(all_results, ignore_index=True)