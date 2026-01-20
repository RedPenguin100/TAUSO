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
        try:
            filenames.remove('Smiles.csv')
        except Exception:
            pass

    if verbose:
        print(f"Loading features from: {filenames}")

    dfs = [pd.read_csv(os.path.join(feature_dir, f)) for f in filenames]

    merged_df = dfs[0]
    for df in dfs[1:]:
        merged_df = pd.merge(merged_df, df, on='index', how='outer')

    return merged_df

def save_feature(df, feature_name):
    feature_dir = SAVED_FEATURES
    os.makedirs(feature_dir, exist_ok=True)
    sub_df = df[['index', feature_name]]
    file_path = os.path.join(feature_dir, f'{feature_name}.csv')
    sub_df.to_csv(file_path, index=False)

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
