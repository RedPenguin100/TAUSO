import pandas as pd
import numpy as np
from notebooks.consts import *
from ...features.vienna_fold import calculate_avg_mfe_over_sense_region


def calculate_mfe_features(df, flank_size=120, window_size=45, step=7):
    """
    Calculates avg MFE for the sense region (+ flanks) based on pre-mRNA coordinates.
    """
    print(f"Starting MFE population with Flank={flank_size}, Window={window_size}, Step={step}")

    df_work = df[[PRE_MRNA_SEQUENCE, SENSE_START, SENSE_LENGTH, 'index']].copy()
    results = []

    for idx, row in df_work.iterrows():

        full_mrna = row[PRE_MRNA_SEQUENCE]
        global_start = row[SENSE_START]
        sense_len = row[SENSE_LENGTH]

        if global_start == -1 or not isinstance(full_mrna, str):
            results.append(np.nan)
            continue

        cut_start = max(0, global_start - flank_size)
        cut_end = min(len(full_mrna), global_start + sense_len + flank_size)

        sub_sequence = full_mrna[cut_start: cut_end]
        relative_sense_start = global_start - cut_start

        if len(sub_sequence) < window_size:
            results.append(np.nan)
            continue

        mfe_score = calculate_avg_mfe_over_sense_region(
            sequence=sub_sequence,
            sense_start=relative_sense_start,
            sense_length=sense_len,
            flank_size=flank_size,
            window_size=window_size,
            step=step
        )

        results.append(mfe_score)

    feature_name = f'mfe_win{window_size}_flank{flank_size}_step{step}'
    result_df = pd.DataFrame({
        'index': df['index'],
        feature_name: results
    })

    return result_df

