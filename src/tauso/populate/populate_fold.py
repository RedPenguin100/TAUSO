import logging
import uuid
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

from .._raccess.core import find_raccess
from ..data.consts import CANONICAL_GENE, SENSE_LENGTH, SENSE_START
from ..features.fold.vienna_fold import calculate_avg_mfe_per_step
from ..features.rna_access.access_calculator import (
    AccessCalculator,
    get_cache,
    get_sense_with_flanks,
)
from ..features.rna_access.rna_access import RNAAccess
from ..features.rna_access.sense_accessibility import compute_sense_accessibility_value
from ..parallel_utils import make_apply_fn

logger = logging.getLogger(__name__)

PRE_MRNA_SEQUENCE = "pre_mrna_sequence"


def validate_cols_in_df(df, cols):
    missing_cols = [c for c in cols if c not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing columns in DataFrame: {missing_cols}")


DEFAULT_SETTINGS = [
    (120, 45, 7),
    (120, 45, 15),
    (120, 45, 10),
    (120, 45, 4),
    (120, 35, 15),
    (120, 35, 10),
    (120, 35, 7),
    (120, 35, 4),
    (120, 55, 15),
    (120, 55, 10),
    (120, 55, 7),
    (120, 55, 4),
    (120, 65, 15),
    (120, 65, 10),
    (120, 65, 7),
    (120, 65, 4),
]


def populate_mfe_features(df, gene_to_data, n_jobs=1, verbose=False, settings=None):
    if settings is None:
        settings = DEFAULT_SETTINGS

    required_cols = [CANONICAL_GENE, SENSE_START, SENSE_LENGTH]
    validate_cols_in_df(df, required_cols)

    unique_genes = df[CANONICAL_GENE].dropna().unique()
    lightweight_gene_to_data = {
        gene: str(gene_to_data[gene].full_mrna) for gene in unique_genes if gene in gene_to_data
    }

    # Group settings by (flank, window). Every entry in a group reuses the same
    # sub-sequence cut and the same set of folded windows; only the `step` grid
    # differs. Doing the apply once across all settings (rather than once per
    # setting) lets us share that work and pays the parallel-dispatch overhead
    # exactly once.
    steps_by_window = defaultdict(list)
    for flank_size, window_size, step in settings:
        steps_by_window[(flank_size, window_size)].append(step)

    feature_names = [f"mfe_win{w}_flank{f}_step{s}" for f, w, s in settings]

    def _process_row(row):
        gene_name = row[CANONICAL_GENE]
        global_start = row[SENSE_START]
        sense_len = row[SENSE_LENGTH]

        out = {name: np.nan for name in feature_names}
        if gene_name not in lightweight_gene_to_data or global_start == -1:
            return out
        full_mrna = lightweight_gene_to_data[gene_name]

        for (flank_size, window_size), steps in steps_by_window.items():
            cut_start = max(0, global_start - flank_size)
            cut_end = min(len(full_mrna), global_start + sense_len + flank_size)
            sub_sequence = full_mrna[cut_start:cut_end]
            if len(sub_sequence) < window_size:
                continue
            step_to_value = calculate_avg_mfe_per_step(
                sequence=sub_sequence,
                sense_start_in_flank=(global_start - cut_start),
                sense_length=sense_len,
                window_size=window_size,
                steps=steps,
            )
            for step, value in step_to_value.items():
                out[f"mfe_win{window_size}_flank{flank_size}_step{step}"] = value
        return out

    apply_fn = make_apply_fn(df, n_jobs=n_jobs, progress_bar=verbose, verbose=2 if verbose else 0)
    results = apply_fn(_process_row, axis=1)
    results_df = pd.DataFrame(list(results), index=df.index)
    for name in feature_names:
        df[name] = results_df[name]
    return df, feature_names


SENSE_AVG_ACCESSIBILITY = "sense_avg_accessibility"

# --- CONSTANTS ---
FLANK_SIZE = 120
ACCESS_SIZE = 13
SEED_SIZE = 13
SEED_SIZES = [SEED_SIZE * m for m in range(1, 4)]
ACCESS_WIN_SIZE = 80


def populate_sense_accessibility(aso_dataframe, gene_to_data):
    access_cache = get_cache(SEED_SIZES, access_size=ACCESS_SIZE)

    valid_rows = [(idx, row) for idx, row in aso_dataframe.iterrows() if row[SENSE_START] != -1]

    for idx, row in valid_rows:
        sense_start = row[SENSE_START]
        sense_length = row[SENSE_LENGTH]
        full_mrna_seq = gene_to_data[row[CANONICAL_GENE]].full_mrna

        flanked_sense = get_sense_with_flanks(full_mrna_seq, sense_start, sense_length, flank_size=FLANK_SIZE)

        if not flanked_sense or len(flanked_sense) < ACCESS_SIZE:
            aso_dataframe.loc[idx, SENSE_AVG_ACCESSIBILITY] = 0
            continue

        avg_sense_access = compute_sense_accessibility_value(
            sense_start,
            sense_length,
            flank=flanked_sense,
            flank_size=FLANK_SIZE,
            access_win_size=ACCESS_WIN_SIZE,
            seed_sizes=SEED_SIZES,
            access_size=ACCESS_SIZE,
            cache=access_cache,
        )
        aso_dataframe.loc[idx, SENSE_AVG_ACCESSIBILITY] = avg_sense_access


DEFAULT_SENSE_CONFIGURATION = [
    {"flank": 120, "access": 20, "seeds": [13]},
    {"flank": 120, "access": 13, "seeds": [4, 6, 8]},
    {"flank": 120, "access": 20, "seeds": [4, 6, 8]},
    {"flank": 120, "access": 13, "seeds": [13, 26, 39]},
]


# =============================================================================
# Sense-accessibility population
# =============================================================================
#
# Goal
# ----
# For each ASO row, compute the mean RNA accessibility of the antisense's target
# window on the pre-mRNA. Output is one feature column per
# (flank_size, access_size, seed_sizes) configuration. Many such configs can be
# evaluated cheaply via `populate_sense_accessibility_multi`.
#
# Pipeline (same for both single-config and multi-config flows)
# -------------------------------------------------------------
# 1. Drop rows with SENSE_START == -1; build a `gene -> full_mrna_string` dict
#    (`_compute_lightweight_gene_to_data`) so worker threads don't pickle/copy
#    the heavyweight gene_to_data objects.
# 2. Split rows into batches; process in parallel threads when n_jobs > 1
#    (raccess is a subprocess, so the GIL doesn't block the workers).
# 3. Per batch:
#    a. `_build_batch_rna_seqs`: slice a flank-padded sense window from each
#       row's mRNA. Returns
#          rna_seqs:        list of (rna_id, RNA-formatted flanked sequence)
#          sense_len_map:   rna_id -> sense window length
#          sense_start_map: rna_id -> sense start inside the flanked window
#       rna_id is the hidden positional `_temp_id`, NOT the original df index;
#       rows whose flanked window is shorter than `min_seq_len` are dropped.
#    b. `RNAAccess.calculate(rna_seqs, seed_sizes, max_span)`: ONE raccess
#       subprocess invocation that emits per-position accessibility energies for
#       every (rna_id, position, seed_size). `max_span` is the per-position
#       folding window raccess uses internally — invariant across configs.
#    c. `_per_position_avg_access` (per config): pick the requested seed_size
#       columns out of the raccess result, slide an `access_size`-wide window
#       (delegated to `calc_access_energies_from_access_res`), and average the
#       per-seed avg columns into a single `avg_access` column.
#    d. `_sense_region_mean`: collapse the per-position frame to one feature
#       value per rna_id by averaging only positions inside the sense window
#       [rel_start, rel_start + sense_length).
# 4. `_assign_feature_column`: write per-batch results back into df_out via the
#    `_temp_id -> value` map. Rows not present in any batch result get NaN.
#
# Why two top-level functions?
# ----------------------------
# `populate_sense_accessibility_batch` handles ONE config and is the historical
# entry point (calculator.py loops over configs and calls it).
# `populate_sense_accessibility_multi` takes a LIST of configs and does (a)-(d)
# in one batched pass: one raccess call per batch with `seed_sizes = sorted
# union of all configs' seeds`, then (c)-(d) once per config. The union call
# costs roughly the same as the most-expensive single-config call (raccess wall
# time is dominated by the smallest seed_size — adding {13,26,39} on top of
# {4,6,8} measured as ~zero cost), so for the 4-config DEFAULT setup this is
# ~3x faster than calling _batch four times.
#
# Both functions go through the EXACT SAME helpers; the diff between them is
# only the dispatch shape (1 config / 1 raccess call vs N configs / 1 union
# call). Output for each config is bit-identical regardless of which function
# you call, because `calc_access_energies_from_access_res` only ever reads its
# requested seed_size's own column in the raccess result.
# =============================================================================


def _build_feature_name(flank_size, access_size, seeds):
    """Stable column name for an accessibility feature. Both _batch and _multi
    call this so the strings stay identical between paths."""
    seeds_str = "-".join(map(str, seeds))
    return f"access_{flank_size}flank_{access_size}access_{seeds_str}seed_sizes"


def _compute_lightweight_gene_to_data(valid_df, gene_to_data):
    """Slim gene_to_data down to a `gene -> mRNA string` dict. Worker threads
    pickle the closure they see; without this step the entire (heavy)
    gene_to_data object gets serialised into every batch."""
    unique_genes = valid_df[CANONICAL_GENE].unique()
    return {gene: str(gene_to_data[gene].full_mrna) for gene in unique_genes if gene in gene_to_data}


def _build_batch_rna_seqs(batch_rows, lightweight_gene_to_data, flank_size, min_seq_len):
    """Cut a `flank_size`-padded window around each row's sense start; format as RNA.

    Rows with a missing gene, or a flanked window shorter than `min_seq_len`
    (typically the smallest access_size in play), are dropped. Returns
    (rna_seqs, sense_len_map, sense_start_map) as described in the section
    header above.
    """
    rna_seqs = []
    sense_len_map = {}
    sense_start_map = {}

    for _, row in batch_rows.iterrows():
        current_id = str(row["_temp_id"])
        gene_name = row[CANONICAL_GENE]
        if gene_name not in lightweight_gene_to_data:
            continue
        full_mrna_seq = lightweight_gene_to_data[gene_name]

        current_sense_start = row[SENSE_START]
        flank_start = max(0, current_sense_start - flank_size)
        relative_start = current_sense_start - flank_start

        flanked_sense = get_sense_with_flanks(
            full_mrna_seq,
            current_sense_start,
            row[SENSE_LENGTH],
            flank_size=flank_size,
        )
        if not flanked_sense or len(flanked_sense) < min_seq_len:
            continue

        flanked_sense = flanked_sense.upper().replace("T", "U")
        rna_seqs.append((current_id, flanked_sense))
        sense_len_map[current_id] = row[SENSE_LENGTH]
        sense_start_map[current_id] = relative_start

    return rna_seqs, sense_len_map, sense_start_map


def _per_position_avg_access(raccess_res, rna_seqs, access_size, seed_sizes):
    """Per-position avg_access frame for ONE (access_size, seed_sizes) config.

    `raccess_res` is the dict returned by `RNAAccess.calculate`: one DataFrame
    of per-position accessibility energies per rna_id, with one column per
    seed_size requested in the raccess call. This helper picks the columns for
    `seed_sizes`, slides an `access_size`-wide window via
    `calc_access_energies_from_access_res` (see that function for the windowing
    math), and averages the per-seed avg columns into a single `avg_access`
    column.

    Returns a DataFrame with columns ["rna_id", "avg_access"] in per-rna_id
    position order. rna_seqs whose length is below `access_size`, and configs
    where every seed > access_size, contribute zero rows.
    """
    eligible_seeds = [s for s in seed_sizes if s <= access_size]
    if not eligible_seeds:
        return pd.DataFrame(columns=["rna_id", "avg_access"])

    df_list = []
    for rna_id, rna_seq in rna_seqs:
        rna_size = len(rna_seq)
        if rna_size < access_size:
            continue
        access_res = raccess_res[rna_id]
        df = AccessCalculator.calc_access_energies_from_access_res(
            access_res,
            rna_size,
            access_size,
            eligible_seeds,
        )
        target_cols = [f"{s}_avg" for s in eligible_seeds if f"{s}_avg" in df.columns]
        df["avg_access"] = df[target_cols].mean(axis=1) if target_cols else 0.0
        df["rna_id"] = rna_id
        df_list.append(df[["rna_id", "avg_access"]])

    if not df_list:
        return pd.DataFrame(columns=["rna_id", "avg_access"])
    return pd.concat(df_list, ignore_index=True)


def _sense_region_mean(per_position_df, sense_len_map, sense_start_map, feature_name):
    """Collapse a per-position avg_access frame to one number per rna_id.

    Averages only positions inside the sense window
    [rel_start, rel_start + sense_length) for each rna_id. Returns DataFrame
    with columns ["_temp_id", <feature_name>].
    """
    if per_position_df.empty:
        return pd.DataFrame(columns=["_temp_id", feature_name])

    df = per_position_df  # caller hands us an owned frame; mutate in place
    df["pos_in_flank"] = df.groupby("rna_id", sort=False, observed=True).cumcount()
    df["sense_length"] = df["rna_id"].map(sense_len_map)
    df["rel_start"] = df["rna_id"].map(sense_start_map)
    mask = (df["pos_in_flank"] >= df["rel_start"]) & (df["pos_in_flank"] < df["rel_start"] + df["sense_length"])

    result = df[mask].groupby("rna_id", as_index=False, observed=True)["avg_access"].mean()
    result.rename(columns={"rna_id": "_temp_id", "avg_access": feature_name}, inplace=True)
    result["_temp_id"] = result["_temp_id"].astype(int)
    return result


def _assign_feature_column(df_out, batch_results, feature_name):
    """Write per-batch results into df_out[`feature_name`] via the _temp_id map.

    `batch_results` is a list of DataFrames with columns ["_temp_id", feature_name]
    (some may be empty). Rows not present in any batch result get NaN.
    """
    non_empty = [br for br in batch_results if not br.empty]
    if not non_empty:
        df_out[feature_name] = np.nan
        return
    results_df = pd.concat(non_empty, ignore_index=True)
    result_map = dict(zip(results_df["_temp_id"], results_df[feature_name]))
    df_out[feature_name] = df_out["_temp_id"].map(result_map)


def populate_sense_accessibility_batch(
    aso_dataframe,
    gene_to_data,
    batch_size=1000,
    flank_size=FLANK_SIZE,
    access_size=ACCESS_SIZE,
    seed_sizes=None,
    access_win_size=ACCESS_WIN_SIZE,
    n_jobs=1,
):
    """Populate ONE accessibility feature column. See section header above for the
    pipeline; this is the single-config form."""
    if seed_sizes is None:
        seed_sizes = SEED_SIZES
    seeds_list = seed_sizes if isinstance(seed_sizes, (list, tuple)) else [seed_sizes]
    feature_name = _build_feature_name(flank_size, access_size, seeds_list)

    df_out = aso_dataframe.copy()
    df_out["_temp_id"] = range(len(df_out))

    valid_mask = df_out[SENSE_START] != -1
    if not valid_mask.any():
        df_out[feature_name] = np.nan
        df_out.drop(columns=["_temp_id"], inplace=True)
        return df_out, feature_name

    valid_df = df_out.loc[valid_mask, ["_temp_id", SENSE_START, SENSE_LENGTH, CANONICAL_GENE]].copy()
    lightweight_gene_to_data = _compute_lightweight_gene_to_data(valid_df, gene_to_data)
    raccess_exe = find_raccess()

    def _process_batch(batch_rows):
        rna_seqs, sense_len_map, sense_start_map = _build_batch_rna_seqs(
            batch_rows,
            lightweight_gene_to_data,
            flank_size,
            min_seq_len=access_size,
        )
        if not rna_seqs:
            return pd.DataFrame(columns=["_temp_id", feature_name])

        # One raccess call for this batch with this config's seeds.
        ra = RNAAccess(seeds_list, access_win_size, raccess_exe)
        ra.set_uuid_for_web(str(uuid.uuid4()))
        raccess_res = ra.calculate(rna_seqs)

        per_position = _per_position_avg_access(raccess_res, rna_seqs, access_size, seeds_list)
        return _sense_region_mean(per_position, sense_len_map, sense_start_map, feature_name)

    batches = [valid_df.iloc[i : i + batch_size] for i in range(0, len(valid_df), batch_size)]
    if n_jobs > 1:
        with ThreadPoolExecutor(max_workers=min(n_jobs, len(batches))) as pool:
            batch_results = list(pool.map(_process_batch, batches))
    else:
        batch_results = [_process_batch(b) for b in batches]

    _assign_feature_column(df_out, batch_results, feature_name)
    df_out.drop(columns=["_temp_id"], inplace=True)
    return df_out, feature_name


def populate_sense_accessibility_multi(
    aso_dataframe,
    gene_to_data,
    configs,
    batch_size=1000,
    access_win_size=ACCESS_WIN_SIZE,
    n_jobs=1,
):
    """Populate ONE feature column per config in `configs` via a single shared
    raccess call per batch. See section header above for the pipeline; this is
    the multi-config form.

    All configs must share `flank` (the per-batch rna_seqs are computed once and
    reused). Output for each config is bit-identical to calling
    populate_sense_accessibility_batch once per config.

    Returns (df, feature_names) with feature_names aligned to `configs`.
    """
    if not configs:
        return aso_dataframe.copy(), []

    flank_sizes = {c["flank"] for c in configs}
    if len(flank_sizes) != 1:
        raise ValueError(f"All configs must share the same flank_size; got {sorted(flank_sizes)}")
    flank_size = next(iter(flank_sizes))

    union_seeds = sorted({s for c in configs for s in c["seeds"]})
    min_access = min(c["access"] for c in configs)
    feature_names = [_build_feature_name(c["flank"], c["access"], c["seeds"]) for c in configs]

    df_out = aso_dataframe.copy()
    df_out["_temp_id"] = range(len(df_out))

    valid_mask = df_out[SENSE_START] != -1
    if not valid_mask.any():
        for fn in feature_names:
            df_out[fn] = np.nan
        df_out.drop(columns=["_temp_id"], inplace=True)
        return df_out, feature_names

    valid_df = df_out.loc[valid_mask, ["_temp_id", SENSE_START, SENSE_LENGTH, CANONICAL_GENE]].copy()
    lightweight_gene_to_data = _compute_lightweight_gene_to_data(valid_df, gene_to_data)
    raccess_exe = find_raccess()

    def _process_batch(batch_rows):
        rna_seqs, sense_len_map, sense_start_map = _build_batch_rna_seqs(
            batch_rows,
            lightweight_gene_to_data,
            flank_size,
            min_seq_len=min_access,
        )
        empty = {fn: pd.DataFrame(columns=["_temp_id", fn]) for fn in feature_names}
        if not rna_seqs:
            return empty

        # ONE raccess call for this batch with the union of all configs' seeds.
        ra = RNAAccess(union_seeds, access_win_size, raccess_exe)
        ra.set_uuid_for_web(str(uuid.uuid4()))
        raccess_res = ra.calculate(rna_seqs)

        return {
            feature_name: _sense_region_mean(
                _per_position_avg_access(raccess_res, rna_seqs, c["access"], c["seeds"]),
                sense_len_map,
                sense_start_map,
                feature_name,
            )
            for c, feature_name in zip(configs, feature_names)
        }

    batches = [valid_df.iloc[i : i + batch_size] for i in range(0, len(valid_df), batch_size)]
    if n_jobs > 1:
        with ThreadPoolExecutor(max_workers=min(n_jobs, len(batches))) as pool:
            batch_results = list(pool.map(_process_batch, batches))
    else:
        batch_results = [_process_batch(b) for b in batches]

    for fn in feature_names:
        _assign_feature_column(df_out, [br[fn] for br in batch_results], fn)
    df_out.drop(columns=["_temp_id"], inplace=True)
    return df_out, feature_names
