import os

import pandas as pd
from joblib import Parallel, delayed

from notebooks.consts import SAVED_FEATURES
from tauso.features.feature_extraction import save_feature_internal


def _get_saved_features_dir(version):
    if version is None:
        return SAVED_FEATURES
    elif version in ["v2", "oligo"]:
        return SAVED_FEATURES.with_name(f"{SAVED_FEATURES.name}_{version}")
    else:
        raise ValueError(f"Unknown version '{version}'")


COMPETITION = [
    "PFRED_PLS",
    "PFRED_SVM",
    "OW_Overall",
    "OW_Tm",
    "OW_Intra_Oligo",
    "OW_Duplex",
    "sfold_accessibility",
    "miranda_score",
    "miranda_energy",
    "oligo_ai_score",
]


def load_all_features(filenames=None, light=True, verbose=False, version=None, n_jobs=1, load_competition=False):
    feature_dir = _get_saved_features_dir(version)
    if not filenames:
        filenames = [f for f in os.listdir(feature_dir) if f.endswith(".csv") and not f.startswith(".")]
        filenames.sort()

    if not filenames:
        raise FileNotFoundError(f"No CSV files found in {feature_dir}")

    if light:
        if "Smiles.csv" in filenames:
            filenames.remove("Smiles.csv")

    if not load_competition:
        filenames = [f for f in filenames if f.removesuffix(".csv") not in COMPETITION]

    if verbose:
        print(f"Loading features from: {filenames}")

    if version is None:
        index = "index"
    else:
        index = "index_" + version

    if n_jobs == 1:
        dfs = [pd.read_csv(os.path.join(feature_dir, f), index_col=index) for f in filenames]
    else:
        dfs = Parallel(n_jobs=n_jobs, backend="threading")(
            delayed(pd.read_csv)(os.path.join(feature_dir, f), index_col=index) for f in filenames
        )

    merged_df = pd.concat(dfs, axis=1, join="outer").reset_index()

    return merged_df


def save_feature(df, feature_name, overwrite=False, version=None):
    save_feature_internal(
        df, feature_name, overwrite=overwrite, version=version, saved_dir_func=_get_saved_features_dir
    )


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
