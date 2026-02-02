import numpy as np
from pandarallel import pandarallel

from ..data.consts import CANONICAL_GENE
from ..features.fold.vienna_fold import calculate_avg_mfe_over_sense_region

SENSE_START = 'sense_start'
SENSE_LENGTH = 'sense_length'
PRE_MRNA_SEQUENCE = 'pre_mrna_sequence'


def validate_cols_in_df(df, cols):
    missing_cols = [c for c in cols if c not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing columns in DataFrame: {missing_cols}")


def populate_mfe_features(df, gene_to_data, n_jobs=1):  # 2. Add n_jobs param
    required_cols = [CANONICAL_GENE, SENSE_START, SENSE_LENGTH]
    validate_cols_in_df(df, required_cols)

    settings = [
        (120, 45, 7), (120, 45, 15), (120, 45, 10), (120, 45, 4),
        (120, 35, 15), (120, 35, 10), (120, 35, 7), (120, 35, 4),
        (120, 55, 15), (120, 55, 10), (120, 55, 7), (120, 55, 4),
        (120, 65, 15), (120, 65, 10), (120, 65, 7), (120, 65, 4),
    ]
    feature_names = []
    for setting in settings:
        flank_size, window_size, step = setting
        print(f"Starting MFE population with Flank={flank_size}, Window={window_size}, Step={step}")

        feature_name = f'mfe_win{window_size}_flank{flank_size}_step{step}'
        feature_names.append(feature_name)

        # TODO: extract outside
        def _process_row(row):
            gene_name = row[CANONICAL_GENE]
            global_start = row[SENSE_START]
            sense_len = row[SENSE_LENGTH]

            if gene_name not in gene_to_data or global_start == -1:
                return np.nan

            locus_info = gene_to_data[gene_name]
            full_mrna = str(locus_info.full_mrna)

            cut_start = max(0, global_start - flank_size)
            cut_end = min(len(full_mrna), global_start + sense_len + flank_size)
            sub_sequence = full_mrna[cut_start: cut_end]

            if len(sub_sequence) < window_size:
                return np.nan

            # Return the score directly
            return calculate_avg_mfe_over_sense_region(
                sequence=sub_sequence,
                sense_start_in_flank=(global_start - cut_start),
                sense_length=sense_len,
                window_size=window_size,
                step=step
            )

        if n_jobs > 1:
            pandarallel.initialize(nb_workers=n_jobs, verbose=0)
            df[feature_name] = df.parallel_apply(_process_row, axis=1)
        else:
            df[feature_name] = df.apply(_process_row, axis=1)

    return df, feature_names
