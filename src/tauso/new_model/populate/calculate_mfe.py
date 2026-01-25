import pandas as pd
import numpy as np
import tauso
_GLOBAL_LOCUS_MAP = None
SENSE_START = 'sense_start'
SENSE_LENGTH = 'sense_length'
PRE_MRNA_SEQUENCE = 'pre_mrna_sequence'


def calculate_mfe_features(df, gene_to_data_map = None, gene_col='canonical_gene_name', flank_size=120, window_size=45, step=7):
    """
    Calculates avg MFE for the sense region (+ flanks) based on pre-mRNA coordinates.
    """
    global _GLOBAL_LOCUS_MAP
    if gene_to_data_map is None:
        if _GLOBAL_LOCUS_MAP is None:
            _GLOBAL_LOCUS_MAP = get_locus_to_data_dict(include_introns=True, genome='GRCh38')

        gene_to_data_map = _GLOBAL_LOCUS_MAP


    print(f"Starting MFE population with Flank={flank_size}, Window={window_size}, Step={step}")

    required_cols = [gene_col, SENSE_START, SENSE_LENGTH, 'index']

    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing columns in DataFrame: {missing_cols}")

    results = []

    # iterate over the dataset rows
    for idx, row in df[required_cols].iterrows():

        gene_name = row[gene_col]
        global_start = row[SENSE_START]
        sense_len = row[SENSE_LENGTH]

        if gene_name not in gene_to_data_map:
            results.append(np.nan)
            continue

        locus_info = gene_to_data_map[gene_name]

        if not hasattr(locus_info, 'full_mrna') or not locus_info.full_mrna:
            results.append(np.nan)
            continue

        full_mrna = str(locus_info.full_mrna)

        if global_start == -1:
            results.append(np.nan)
            continue

        # Slicing the relevant sequence out of the RNA
        cut_start = max(0, global_start - flank_size)
        cut_end = min(len(full_mrna), global_start + sense_len + flank_size)

        sub_sequence = full_mrna[cut_start: cut_end]
        relative_sense_start = global_start - cut_start

        if len(sub_sequence) < window_size:
            results.append(np.nan)
            continue

        # MFE calculation
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

