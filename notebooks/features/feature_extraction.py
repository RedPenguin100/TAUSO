import pandas as pd
import os
from notebooks.consts import SAVED_FEATURES


def load_all_features(filenames=None, light=True, verbose=False):
    # this function loads all the features and the current data, merges all the features according to the index
    # and returns a merged df with all the features and the data
    # missing values would be considered as NaN!

    feature_dir = SAVED_FEATURES
    if not filenames:
        filenames = [f for f in os.listdir(feature_dir) if f.endswith('.csv') and not f.startswith('.')]
        filenames.sort()

    if not filenames:
        raise FileNotFoundError(f"No CSV files found in {feature_dir}")

    if light:
        filenames.remove('Smiles.csv')

    if verbose:
        print(f"Loading features from: {filenames}")

    dfs = [pd.read_csv(os.path.join(feature_dir, f)) for f in filenames]

    merged_df = dfs[0]
    for df in dfs[1:]:
        merged_df = pd.merge(merged_df, df, on='index', how='outer')

    return merged_df


import os
import pandas as pd
import numpy as np


def _check_conflicts(new_sub_df, file_path, feature_name):
    existing_df = pd.read_csv(file_path)

    # Merge on 'index' to compare aligned rows
    comparison = pd.merge(
        new_sub_df,
        existing_df,
        on='index',
        suffixes=('_new', '_existing'),
        how='inner'
    )

    col_new = f'{feature_name}_new'
    col_old = f'{feature_name}_existing'

    # --- FIX START: Handle Floating Point Precision ---
    # Check if the column is numeric (float/int)
    if pd.api.types.is_numeric_dtype(comparison[col_new]):
        # Option A: Round both to 6 decimals (easiest to read)
        # Option B: Use np.isclose for math precision

        # We create a mask where values are NOT close
        is_close = np.isclose(comparison[col_new], comparison[col_old], rtol=1e-05, atol=1e-08, equal_nan=True)
        diff_mask = ~is_close  # Invert it: True where NOT close
    else:
        # For strings/other types, keep strict comparison
        diff_mask = comparison[col_new] != comparison[col_old]
    # --- FIX END ---

    diff_rows = comparison[diff_mask]

    if not diff_rows.empty:
        print(f"!!! Conflict detected for feature: {feature_name} !!!")
        print(f"Found {len(diff_rows)} differing rows. Top 10 differences:")
        print(diff_rows[['index', col_new, col_old]].head(10))
        print("-" * 30)
        raise FileExistsError(
            f"Feature '{feature_name}' exists and has conflicting values. "
            "Pass overwrite=True to force save."
        )

    print(f"File exists for '{feature_name}' but values are identical (within tolerance). No action taken.")


# The main function remains the same
def save_feature(df, feature_name, overwrite=False):
    feature_dir = SAVED_FEATURES
    os.makedirs(feature_dir, exist_ok=True)

    sub_df = df[['index', feature_name]].copy()
    file_path = os.path.join(feature_dir, f'{feature_name}.csv')

    file_exists = os.path.exists(file_path)

    if not file_exists or overwrite:
        sub_df.to_csv(file_path, index=False)
        action = "Overwrote" if file_exists else "Saved"
        print(f"{action} feature: {feature_name}")
    else:
        _check_conflicts(sub_df, file_path, feature_name)

# TODO: fix this function
# def read_base_df():
#     all_data = pd.read_csv(DATA_PATH_NEW / 'data_asoptimizer_updated.11.7.csv',
#                            dtype={'index': int, 'ISIS': int, 'Target_gene': str,
#                                   'Cell_line': str, 'Density(cells_per_well)': float,
#                                   'Transfection': str, 'ASO_volume(nm)': float, 'Treatment_Period(hours)': float,
#                                   'Primer_probe_set': str, 'Sequence': str, 'Modification': str, 'Location': str,
#                                   'Chemical_Pattern': str, 'Linkage': str, 'Linage_Location' : str, 'Smiles' : str,
#                                   'Inhibition(%)' : float, 'seq_length' : int, 'Canonical Gene Name' : str,
#                                   'Cell line organism' : str, 'Transcript' : str, 'Location_in_sequence' : int,
#                                   'Location_div_by_length' : float, 'true_length_of_seq' : int, 'mod_scan' : int,
#                                   'cell_line_uniform' : str
#                                   })
#     return all_data
