import pyBigWig
import numpy as np

from ..features.context.ribo_seq import RIBOSEQ_40S_HUMAN_DATA, calculate_ribo_seq_row


def populate_ribo_seq(organism, aso_df, flanks=(0, 10, 20, 50, 100, 125, 150), how="mean"):
    if organism != 'human':
        raise ValueError("Unsupported organism for ribo_seq feature")

    if how not in {"sum", "mean", "max", "nz_mean"}:
        raise ValueError(how)

    bw = pyBigWig.open(str(RIBOSEQ_40S_HUMAN_DATA))

    # Pre-initialize columns with NaN to ensure structure exists even if rows fail
    # We define the column names first based on flanks and 'how'
    new_features_list = [f"ribo_gene_{how}"]
    for f in flanks:
        new_features_list.append(f"ribo_f{f}_{how}")
        if f > 0:
            new_features_list.append(f"ribo_upstream_{f}_{how}")
            new_features_list.append(f"ribo_downstream_{f}_{how}")

    # Initialize new columns in the original dataframe
    for col in new_features_list:
        aso_df[col] = np.nan

    # Iterate by index to assign directly in place
    for idx, row in aso_df.iterrows():
        try:
            feat_dict = calculate_ribo_seq_row(row, bw, flanks, how)

            # Assign values directly to the specific index
            for key, val in feat_dict.items():
                aso_df.at[idx, key] = val

        except Exception as e:
            print(f"Row {idx} failed: {e}")
            continue

    return aso_df, new_features_list
