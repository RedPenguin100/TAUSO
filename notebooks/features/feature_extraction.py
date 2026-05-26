import concurrent.futures
import os

import pandas as pd

from notebooks.consts import SAVED_FEATURES
from tauso.populate.feature_cache import save_feature_internal


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
    "ClinASO_score",
]


def get_dtype_for_feature(filename, index_col_name):
    if filename.startswith("RBP_"):
        feat_type = "float64"
    elif filename.startswith("OHE_"):
        feat_type = "int8"
    elif filename.startswith("tAI_") or filename.startswith("CAI_") or filename.startswith("ENC_"):
        feat_type = "float64"
    elif filename.startswith("access_") or filename.startswith("mfe"):
        feat_type = "float64"
    elif filename.startswith("off_target_score_") or filename.startswith("off_target_single_"):
        feat_type = "float64"
    elif filename in [
        "RNaseH1_Inefficacy_GGGG.csv",
        "RNaseH1_Potency_TCCC.csv",
        "RNaseH1_Potency_TCTC.csv",
        "RNaseH1_Potency_TTCC.csv",
    ]:
        feat_type = "int"
    elif filename.startswith("RNaseH1_Krel") or filename.startswith("RNaseH1_score"):
        feat_type = "float64"

    elif filename.startswith("Modification_"):
        feat_type = "float64"
    elif filename in [
        "MOE_DIFF_37_MD_GB_HYBR.csv",
        "MOE_DIFF_37_MD_PB_HYBR.csv",
        "METHYL_CYTOSINES.csv",
        "LNA_DIFF_37_HYBR.csv",
        "CET_DIFF_37_HYBR.csv",
        "DNA_HYBR_DIFF.csv",
        "PSDNA_RNA_MD_37_GB_TOTAL_HYBR.csv",
        "PSDNA_RNA_MD_37_PB_TOTAL_HYBR.csv",
        "PSDNA_DIFF_37_HYBR.csv",
        "TOTAL_DNA_HYBR.csv",
        "TOTAL_DNA_RNA_HYBR.csv",
        "TOTAL_PSDNA_HYBR.csv",
    ]:
        feat_type = "float64"
    elif filename.startswith("Sequence_"):
        feat_type = "float64"
    elif filename.startswith("SeqChem_"):
        feat_type = "float64"
    elif filename.startswith("ribo_"):
        feat_type = "float64"
    elif filename in [
        "sense_exon.csv",
        "sense_intron.csv",
        "sense_length.csv",
        "sense_start.csv",
        "sense_start_from_end.csv",
        "sense_utr.csv",
        "sense_3utr.csv",
        "sense_5utr.csv",
        "n_sense_exon.csv",
        "n_sense_intron.csv",
        "n_sense_utr.csv",
        "n_sense_3utr.csv",
        "n_sense_5utr.csv",
    ]:
        feat_type = "int"
    elif filename in ["on_target_total_hybridization_0.csv", "on_target_total_hybridization_1200.csv"]:
        feat_type = "float64"
    elif filename in [
        "target_expression.csv",
        "rnase_expression.csv",
        "stab2_expression.csv",
        "mrc1_expression.csv",
        "msr1_expression.csv",
    ]:
        feat_type = "float64"
    elif filename in ["max_consecutive_PO.csv", "ps_end_score.csv"]:
        feat_type = "int"
    elif filename == "po_percentage.csv":
        feat_type = "float64"
    elif filename in [
        "chem_2nd_gen.csv",
        "chem_3rd_gen.csv",
        "Electroporation.csv",
        "Gymnosis.csv",
        "Lipofection.csv",
        "Other.csv",
        "is_hepa.csv",
        "is_all_ps.csv",
    ]:
        feat_type = "int"
    elif filename in ["sense_type.csv"]:
        feat_type = "str"
    elif filename in ["HalfLife_Source.csv", "Mapped_Cell_Proxy.csv"]:
        feat_type = "str"
    elif filename in ["mRNA_HalfLife.csv"]:
        feat_type = "float64"
    # Competition
    elif filename.startswith("PFRED_") or filename.startswith("OW_"):
        feat_type = "float64"
    elif filename.startswith("sfold_") or filename.startswith("oligo_ai_") or filename.startswith("miranda_"):
        feat_type = "float64"
    elif filename.startswith("ClinASO"):
        feat_type = "float64"
    else:
        print(f"Filename not assigned: {filename}")
        feat_type = None
    col_name = filename.removesuffix(".csv")

    if feat_type is not None:
        return {index_col_name: "int32", col_name: feat_type}
    return {index_col_name: "int32"}


def load_all_features(filenames=None, light=True, verbose=False, version=None, load_competition=False):
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

    index_col_name = "index" if version is None else f"index_{version}"

    # 2. Create a tiny helper function to read a single file
    def _read_single_csv(f):
        return pd.read_csv(
            os.path.join(feature_dir, f),
            index_col=index_col_name,
            dtype=get_dtype_for_feature(f, index_col_name),
            engine="pyarrow",
            dtype_backend="pyarrow",
        )

    # 3. Fire off the reads in parallel!
    # This automatically spins up multiple threads to read your CSVs simultaneously.
    # It returns the 'dfs' list in the exact same order as your original for-loop.
    with concurrent.futures.ThreadPoolExecutor() as executor:
        dfs = list(executor.map(_read_single_csv, filenames))

    # 4. Concat remains untouched
    merged_df = pd.concat(dfs, axis=1, join="outer", copy=False).reset_index()

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
