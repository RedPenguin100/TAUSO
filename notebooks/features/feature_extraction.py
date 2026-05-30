import concurrent.futures
import logging
import os

import pandas as pd
import pyarrow.parquet as pq

from tauso.populate.feature_cache import (
    FEATURE_CACHE_FILES,
    ensure_cache,
    feature_store_dir,
    save_feature_internal,
)

logger = logging.getLogger(__name__)


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
    elif filename.startswith("seq_"):
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


def _index_col(version):
    if version in FEATURE_CACHE_FILES:
        return FEATURE_CACHE_FILES[version]["index_col"]
    return f"index_{version}"


def _maybe_fetch_cache(version, feature_dir):
    """If `version` is a registered run, ensure the wide cache is in `feature_dir`; return its path."""
    if version not in FEATURE_CACHE_FILES:
        return None
    cache_file = os.path.join(feature_dir, FEATURE_CACHE_FILES[version]["filename"])
    if not os.path.exists(cache_file):
        ensure_cache(version)
    return cache_file


def _feature_name(path):
    return os.path.basename(path).removesuffix(".parquet").removesuffix(".csv")


def _keep_filter(light, load_competition, wanted):
    excluded = set() if load_competition else set(COMPETITION)
    if light:
        excluded.add("Smiles")
    wanted_set = set(wanted) if wanted else None
    return lambda name: name not in excluded and (wanted_set is None or name in wanted_set)


def _list_loose_files(feature_dir, cache_file):
    """Loose per-feature files (.csv or .parquet) in `feature_dir`, excluding both
    the active wide cache and any *other* wide-cache parquets left over from
    previous versions.

    Active wide cache: filename matches the one registered in FEATURE_CACHE_FILES.
    Stale wide cache: any other .parquet whose column count is wider than a
    single-feature file (heuristic: more than `WIDE_CACHE_COL_THRESHOLD` columns).
    Loading a stale wide cache alongside the active one silently overrides the
    active cache's values for every shared column (loose files take precedence in
    load_all_features), so we refuse and warn instead.
    """
    cache_name = os.path.basename(cache_file) if cache_file else None
    WIDE_CACHE_COL_THRESHOLD = 5  # single-feature files have 2 cols (index + value)
    out = []
    stale_caches = []
    for f in sorted(os.listdir(feature_dir)):
        if f.startswith("."):
            continue
        if not (f.endswith(".csv") or f.endswith(".parquet")):
            continue
        if f == cache_name:
            continue
        full_path = os.path.join(feature_dir, f)
        if f.endswith(".parquet"):
            try:
                n_cols = len(pq.read_schema(full_path).names)
            except Exception:
                n_cols = 0
            if n_cols > WIDE_CACHE_COL_THRESHOLD:
                stale_caches.append((f, n_cols))
                continue
        out.append(full_path)

    if stale_caches:
        details = ", ".join(f"{n} ({c} cols)" for n, c in stale_caches)
        logger.warning(
            "Skipping %d stale wide-cache parquet(s) in %s: %s. "
            "These look like leftover caches from previous FEATURE_CACHE_FILES versions; "
            "loading them alongside the current cache would silently override its values "
            "for any shared column. Delete or move them aside to silence this warning.",
            len(stale_caches), feature_dir, details,
        )
    return out


def _read_loose_file(path, index_col):
    """Read a loose per-feature file (CSV or Parquet) as a DataFrame indexed by `index_col`."""
    if path.endswith(".csv"):
        return pd.read_csv(
            path,
            index_col=index_col,
            dtype=get_dtype_for_feature(os.path.basename(path), index_col),
            engine="pyarrow",
            dtype_backend="pyarrow",
        )
    return pd.read_parquet(path).set_index(index_col)


def _read_cache(cache_file, index_col, excluded_cols, keep_name):
    """Read the wide cache, skipping columns covered by shards or filtered out."""
    cols = [c for c in pq.read_schema(cache_file).names if c != index_col and keep_name(c) and c not in excluded_cols]
    if not cols:
        return None
    return pd.read_parquet(cache_file, columns=[index_col] + cols).set_index(index_col)


def load_all_features(filenames=None, light=True, verbose=False, version="oligo", load_competition=False):
    feature_dir = _get_saved_features_dir(version)
    cache_file = _maybe_fetch_cache(version, feature_dir)

    keep = _keep_filter(light=light, load_competition=load_competition, wanted=filenames)
    loose_files = [p for p in _list_loose_files(feature_dir, cache_file) if keep(_feature_name(p))]
    idx = _index_col(version)

    with concurrent.futures.ThreadPoolExecutor() as ex:
        frames = list(ex.map(lambda p: _read_loose_file(p, idx), loose_files))

    if cache_file:
        loose_cols = {c for f in frames for c in f.columns}
        cache_frame = _read_cache(cache_file, idx, excluded_cols=loose_cols, keep_name=keep)
        if cache_frame is not None:
            frames.append(cache_frame)

    if not frames:
        raise FileNotFoundError(f"No features found in {feature_dir}")
    merged = pd.concat(frames, axis=1, join="outer", copy=False).reset_index()
    if verbose:
        print(f"Loaded {merged.shape[1] - 1} feature columns from {feature_dir}")
    return merged


def iter_all_features(batch_size=5000, light=True, version="oligo", load_competition=False, filenames=None):
    """Stream the merged feature matrix in row batches via pyarrow (peak ~ one batch).

    Requires a registered cache for `version`; loose per-feature files are sliced per batch.
    """
    feature_dir = _get_saved_features_dir(version)
    cache_file = _maybe_fetch_cache(version, feature_dir)
    if cache_file is None:
        raise FileNotFoundError(f"iter_all_features needs a registered cache for run '{version}'")

    keep = _keep_filter(light=light, load_competition=load_competition, wanted=filenames)
    idx = _index_col(version)
    loose_paths = [p for p in _list_loose_files(feature_dir, cache_file) if keep(_feature_name(p))]
    loose_data = {_feature_name(p): _read_loose_file(p, idx)[_feature_name(p)] for p in loose_paths}

    cols = [c for c in pq.read_schema(cache_file).names if c != idx and keep(c) and c not in loose_data]
    parquet = pq.ParquetFile(cache_file)
    for batch in parquet.iter_batches(batch_size=batch_size, columns=[idx] + cols):
        bdf = batch.to_pandas().set_index(idx)
        for name, col in loose_data.items():
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
