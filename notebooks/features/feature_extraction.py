import concurrent.futures
import logging
import os

import pandas as pd
import pyarrow.parquet as pq

from tauso.populate.feature_cache import (
    FEATURE_CACHE_FILES,
    LOOSE_SHARD_SUBDIR,
    ZENODO_COMPETITION_SHARDS,
    ensure_cache,
    ensure_competition_cache,
    feature_store_dir,
    loose_shard_dir,
    save_feature_internal,
)

logger = logging.getLogger(__name__)
_warned_legacy_shards = set()


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
    "OW_Break_Target",
    "sfold_accessibility",
    "miranda_score",
    "miranda_energy",
    "oligo_ai_score",
    "ClinASO_score",
]


def get_dtype_for_feature(filename, index_col_name):
    name = filename.removesuffix(".parquet").removesuffix(".csv")
    if name.startswith("rbp_"):
        feat_type = "float64"
    elif name.startswith("ohe_"):
        feat_type = "int8"
    elif name.startswith("tai_") or name.startswith("cai_") or name.startswith("enc_"):
        feat_type = "float64"
    elif name.startswith("fold_"):
        feat_type = "float64"
    elif name.startswith("off_target_score_") or name.startswith("off_target_single_"):
        feat_type = "float64"
    elif name.startswith("rnase_"):
        feat_type = "float64"

    elif name in {
        # mod_ps_* counts/booleans are int; the rest of the mod_* family (mod_sugar_* +
        # mod_ps_po_percentage / mod_ps_frac_mod / mod_ps_frac_dna) is float64.
        "mod_ps_max_consecutive_po",
        "mod_ps_end_score",
        "mod_ps_wing5_count",
        "mod_ps_gap_count",
        "mod_ps_wing3_count",
    }:
        feat_type = "int"
    elif name.startswith("mod_"):
        feat_type = "float64"
    elif name.startswith("hybr_"):
        feat_type = "float64"
    elif name.startswith("seq_"):
        feat_type = "float64"
    elif name.startswith("ribo_"):
        feat_type = "float64"
    elif name.startswith("expr_"):
        feat_type = "float64"
    elif name in {
        # tox_g4hunter_max is a max windowed G4Hunter score (real-valued); tox_3prime_g_fraction
        # is a fraction in [0, 1]. Both are float; everything else in the tox_* family is a
        # motif/length count.
        "tox_g4hunter_max",
        "tox_3prime_g_fraction",
    }:
        feat_type = "float64"
    elif name.startswith("term5p_") or name.startswith("tox_"):
        feat_type = "int"
    elif name in {"structure_sense_length", "structure_sense_start", "structure_sense_start_from_end"}:
        feat_type = "int"
    elif name in {
        "struct_sense_in_exon",
        "struct_sense_in_intron",
        "struct_sense_in_utr",
        "struct_sense_in_3utr",
        "struct_sense_in_5utr",
        "struct_sense_in_cds",
        "struct_sense_in_cds_non_exclusive",
        "struct_sense_in_exon_non_exclusive",
    }:
        feat_type = "int8"
    elif name in {"structure_sense_start_norm", "structure_sense_start_from_end_norm"}:
        feat_type = "float64"
    elif name.startswith("structure_sense_dist_to_") or name.startswith("structure_sense_mrna_dist_to_"):
        feat_type = "float64"
    elif name.startswith("on_target_total_hybridization_"):
        feat_type = "float64"
    elif name in {
        "chem_1st_gen",
        "chem_2nd_gen",
        "chem_3rd_gen",
    }:
        feat_type = "int"
    elif name in {
        # Transfection one-hots are float64 because rows with an unrecognized
        # transfection_raw value (literal "Other", missing, or any label outside
        # the three) carry NaN across all three columns -- NaN cannot live in
        # an int column.
        "transfection_electroporation",
        "transfection_gymnosis",
        "transfection_lipofection",
    }:
        feat_type = "float64"
    elif name == "structure_sense_type":
        feat_type = "str"
    elif name in {"halflife_source", "halflife_cell_proxy"}:
        feat_type = "str"
    elif name == "halflife_value":
        feat_type = "float64"
    # Competition
    elif name.startswith("PFRED_") or name.startswith("OW_"):
        feat_type = "float64"
    elif name.startswith("sfold_") or name.startswith("oligo_ai_") or name.startswith("miranda_"):
        feat_type = "float64"
    elif filename.startswith("ClinASO"):
        feat_type = "float64"
    else:
        logger.warning("Filename not assigned: %s", filename)
        feat_type = None

    if feat_type is not None:
        return {index_col_name: "int32", name: feat_type}
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


def _maybe_fetch_competition(version, load_competition):
    """Fetch the competition shards from Zenodo into _patches/ when requested and registered."""
    if load_competition and version in ZENODO_COMPETITION_SHARDS:
        ensure_competition_cache(version)


def _feature_name(path):
    return os.path.basename(path).removesuffix(".parquet").removesuffix(".csv")


def _keep_filter(light, load_competition, wanted):
    excluded = set() if load_competition else set(COMPETITION)
    if light:
        excluded.add("Smiles")
    wanted_set = set(wanted) if wanted else None
    return lambda name: name not in excluded and (wanted_set is None or name in wanted_set)


def _list_loose_files(feature_dir, cache_file):
    """Loose per-feature shards (.parquet/.csv). Canonical location is `<feature_dir>/_patches/`.

    During the migration window we still pick up shards left in the bare feature_dir, and
    skip wide-cache parquets there via a column-count heuristic. A one-shot warning per
    directory tells the user how to move legacy shards into _patches/.
    """
    patches_dir = loose_shard_dir(feature_dir)
    out = []
    if os.path.isdir(patches_dir):
        for f in sorted(os.listdir(patches_dir)):
            if f.startswith(".") or not (f.endswith(".parquet") or f.endswith(".csv")):
                continue
            out.append(os.path.join(patches_dir, f))
    patches_names = {_feature_name(p) for p in out}

    cache_name = os.path.basename(cache_file) if cache_file else None
    WIDE_CACHE_COL_THRESHOLD = 5  # single-feature files have 2 cols (index + value)
    legacy = []
    stale_caches = []
    for f in sorted(os.listdir(feature_dir)):
        if f.startswith("."):
            continue
        if not (f.endswith(".csv") or f.endswith(".parquet")):
            continue
        if f == cache_name:
            continue
        full_path = os.path.join(feature_dir, f)
        if _feature_name(full_path) in patches_names:
            continue  # superseded by the canonical _patches shard
        if f.endswith(".parquet"):
            try:
                n_cols = len(pq.read_schema(full_path).names)
            except Exception:
                n_cols = 0
            if n_cols > WIDE_CACHE_COL_THRESHOLD:
                stale_caches.append((f, n_cols))
                continue
        legacy.append(full_path)

    if legacy and feature_dir not in _warned_legacy_shards:
        _warned_legacy_shards.add(feature_dir)
        logger.warning(
            "Found %d legacy loose shard(s) in %s. New writes go to %s/. "
            "To migrate: mkdir -p %s && mv %s/*.parquet %s/",
            len(legacy), feature_dir, LOOSE_SHARD_SUBDIR, patches_dir, feature_dir, patches_dir,
        )
    out.extend(legacy)

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
    _maybe_fetch_competition(version, load_competition)

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
    _maybe_fetch_competition(version, load_competition)

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
