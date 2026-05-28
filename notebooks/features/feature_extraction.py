import os

import pandas as pd
import pyarrow.parquet as pq

from tauso.populate.feature_cache import (
    cache_path,
    ensure_cache,
    feature_store_dir,
    save_feature_internal,
)


def _get_saved_features_dir(version):
    if version is None:
        raise ValueError("version is required to resolve the feature store directory")
    return feature_store_dir(version)


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


def get_dtype_for_feature(filename, index_col_name):
    if filename.startswith("RBP_"):
        feat_type = "float64"
    elif filename.startswith("OHE_"):
        feat_type = "int8"
    elif filename.startswith("tAI_") or filename.startswith("CAI_") or filename.startswith("ENC_"):
        feat_type = "float64"
    elif filename.startswith("fold_"):
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
    elif filename.startswith("hybr_"):
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
    else:
        print(f"Filename not assigned: {filename}")
        feat_type = None
    col_name = filename.removesuffix(".csv")

    if feat_type is not None:
        return {index_col_name: "int32", col_name: feat_type}
    return {index_col_name: "int32"}


def load_all_features(filenames=None, light=True, verbose=False, version="oligo", load_competition=False):
    feature_dir = feature_store_dir(version)
    index_col = f"index_{version}"

    # The wide frozen cache is just another parquet in the folder; fetch it if it isn't here yet.
    cache_file = cache_path(version)
    if not os.path.exists(cache_file):
        try:
            ensure_cache(version)
        except FileNotFoundError:
            # No published cache and no local copy: fall back to dev shards, if any exist.
            if not any(f.endswith(".parquet") for f in os.listdir(feature_dir)):
                raise

    parquet_files = sorted(os.path.join(feature_dir, f) for f in os.listdir(feature_dir) if f.endswith(".parquet"))
    if not parquet_files:
        raise FileNotFoundError(f"No parquet features found in {feature_dir}")

    excluded = set() if load_competition else set(COMPETITION)
    if light:
        excluded.add("Smiles")
    wanted = set(filenames) if filenames else None

    def keep_name(name):
        return name not in excluded and (wanted is None or name in wanted)

    # Loose per-feature dev shards take precedence over the cache on a name clash.
    frames = []
    overridden = set()
    for path in parquet_files:
        if path == cache_file:
            continue
        name = os.path.basename(path).removesuffix(".parquet")
        if not keep_name(name):
            continue
        frames.append(pd.read_parquet(path).set_index(index_col))
        overridden.add(name)

    # Read the cache (if present), skipping excluded / shard-overridden / unwanted columns.
    if os.path.exists(cache_file):
        cache_cols = [c for c in pq.read_schema(cache_file).names if c != index_col]
        keep = [c for c in cache_cols if keep_name(c) and c not in overridden]
        if keep:
            frames.append(pd.read_parquet(cache_file, columns=[index_col] + keep).set_index(index_col))

    if not frames:
        raise FileNotFoundError(f"No matching features found in {feature_dir}")

    merged_df = pd.concat(frames, axis=1, join="outer", copy=False).reset_index()
    if verbose:
        print(
            f"Loaded {merged_df.shape[1] - 1} feature columns from {feature_dir} ({len(overridden)} dev shard(s) + cache)"
        )
    return merged_df


def iter_all_features(batch_size=5000, light=True, version="oligo", load_competition=False, filenames=None):
    """Stream the merged feature matrix in row batches, for low peak memory on small machines.

    Yields DataFrames of up to ``batch_size`` rows with the same columns as load_all_features.
    The wide cache is read row-group by row-group via pyarrow (peak ~= one batch); the small
    per-feature dev shards are held in full and sliced per batch (and still override the cache).
    Requires the cache to be present (the cache is the large thing worth streaming).
    """
    feature_dir = feature_store_dir(version)
    index_col = f"index_{version}"
    cache_file = cache_path(version)
    if not os.path.exists(cache_file):
        ensure_cache(version)

    excluded = set() if load_competition else set(COMPETITION)
    if light:
        excluded.add("Smiles")
    wanted = set(filenames) if filenames else None

    def keep_name(name):
        return name not in excluded and (wanted is None or name in wanted)

    # Dev shards are small (one column each): load fully, slice per batch; they override the cache.
    shards = {}
    for f in sorted(os.listdir(feature_dir)):
        path = os.path.join(feature_dir, f)
        if f.endswith(".parquet") and path != cache_file:
            name = f.removesuffix(".parquet")
            if keep_name(name):
                shards[name] = pd.read_parquet(path).set_index(index_col)[name]

    cache_cols = [c for c in pq.read_schema(cache_file).names if c != index_col]
    keep = [c for c in cache_cols if keep_name(c) and c not in shards]

    parquet = pq.ParquetFile(cache_file)
    for batch in parquet.iter_batches(batch_size=batch_size, columns=[index_col] + keep):
        bdf = batch.to_pandas().set_index(index_col)
        for name, col in shards.items():
            bdf[name] = col.reindex(bdf.index)
        yield bdf.reset_index()


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
