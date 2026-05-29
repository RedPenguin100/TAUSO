import logging
import uuid
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

from .._raccess.core import find_raccess
from ..data.consts import CANONICAL_GENE, SENSE_LENGTH, SENSE_START
from ..features.fold.vienna_fold import calculate_avg_mfe_per_step
from ..features.rna_access.access_calculator import get_sense_with_flanks, window_access_energies
from ..features.rna_access.rna_access import RNAAccess
from ..parallel_utils import make_apply_fn

logger = logging.getLogger(__name__)


def validate_cols_in_df(df, cols):
    missing_cols = [c for c in cols if c not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing columns in DataFrame: {missing_cols}")


# MFE grid: (flank, window, step). Spans tight-local (flank=30) through
# medium (60), wide (120), and global (240) folding contexts.
DEFAULT_SETTINGS = [
    # Tight local (flank=30)
    (30, 20, 4),
    (30, 30, 7),
    (30, 45, 7),
    (30, 65, 7),
    (30, 65, 15),
    (30, 90, 15),
    # Medium context
    (60, 55, 7),
    # Wide context
    (120, 55, 7),
    (120, 65, 7),
    # Global context — captures long-range stems
    (240, 90, 15),
]


def populate_mfe_features(df, gene_to_data, n_jobs=1, verbose=False, settings=None):
    if settings is None:
        settings = DEFAULT_SETTINGS

    required_cols = [CANONICAL_GENE, SENSE_START, SENSE_LENGTH]
    validate_cols_in_df(df, required_cols)

    lightweight_gene_to_data = _lightweight_gene_to_data(df[CANONICAL_GENE].dropna().unique(), gene_to_data)

    # Group settings by (flank, window). Every entry in a group reuses the same
    # sub-sequence cut and the same set of folded windows; only the `step` grid
    # differs. Doing the apply once across all settings (rather than once per
    # setting) lets us share that work and pays the parallel-dispatch overhead
    # exactly once.
    steps_by_window = defaultdict(list)
    for flank_size, window_size, step in settings:
        steps_by_window[(flank_size, window_size)].append(step)

    feature_names = [f"fold_mfe_win{w}_flank{f}_step{s}" for f, w, s in settings]

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
                out[f"fold_mfe_win{window_size}_flank{flank_size}_step{step}"] = value
        return out

    apply_fn = make_apply_fn(df, n_jobs=n_jobs, progress_bar=verbose, verbose=2 if verbose else 0)
    results = apply_fn(_process_row, axis=1)
    results_df = pd.DataFrame(list(results), index=df.index)
    for name in feature_names:
        df[name] = results_df[name]
    return df, feature_names


# Default config for the single-config entry point (populate_sense_accessibility_batch).
FLANK_SIZE = 120
ACCESS_SIZE = 13
SEED_SIZES = [13, 26, 39]
ACCESS_WIN_SIZE = 80


# raccess grid: (flank, access, seeds). Spans tight-local (flank=15) through
# wide (flank=120). Capped at flank=120 because raccess scales O(L²)/O(L³)
# in window length and flank=240 would be prohibitive.
DEFAULT_SENSE_CONFIGURATION = [
    # Tight local (flank=15) — immediate stem-loop occlusion
    {"flank": 15, "access": 8, "seeds": [4, 6, 8]},
    {"flank": 15, "access": 13, "seeds": [4, 6, 8]},
    # Local (flank=30) — short-range context
    {"flank": 30, "access": 13, "seeds": [4, 6, 8]},
    {"flank": 30, "access": 20, "seeds": [4, 6, 8]},
    # Short regional (flank=60)
    {"flank": 60, "access": 13, "seeds": [4, 6, 8]},
    {"flank": 60, "access": 20, "seeds": [4, 6, 8]},
    # Wider (flank=120)
    {"flank": 120, "access": 6, "seeds": [4, 6]},
    {"flank": 120, "access": 8, "seeds": [4, 6, 8]},
    {"flank": 120, "access": 10, "seeds": [4, 6, 8]},
    {"flank": 120, "access": 13, "seeds": [4, 6, 8]},
    {"flank": 120, "access": 13, "seeds": [13, 26, 39]},
    {"flank": 120, "access": 20, "seeds": [4, 6, 8]},
    {"flank": 120, "access": 39, "seeds": [13, 26, 39]},
]


# =============================================================================
# Sense-accessibility population
# =============================================================================
# For each ASO row we compute the mean RNA accessibility of its target window on
# the pre-mRNA, one feature column per (flank, access_size, seeds) config.
#
# Per batch of rows: cut a flank-padded window around each sense start, run ONE
# raccess call to get per-position opening energies, slide an access_size window
# over them (window_access_energies), then average across the sense region.
#
# populate_sense_accessibility_multi is the core (one raccess call per batch over
# the union of all configs' seeds; raccess wall time is dominated by the smallest
# seed, so extra seeds are nearly free). populate_sense_accessibility_batch is a
# thin single-config wrapper over it.
# =============================================================================


def _build_feature_name(flank_size, access_size, seeds):
    """Stable column name for an accessibility config."""
    return f"fold_access_{flank_size}flank_{access_size}access_{'-'.join(map(str, seeds))}seed_sizes"


def _lightweight_gene_to_data(genes, gene_to_data):
    """gene -> mRNA string, so worker threads don't pickle the heavy gene_to_data."""
    return {gene: str(gene_to_data[gene].full_mrna) for gene in genes if gene in gene_to_data}


def _build_batch_rna_seqs(batch_rows, lightweight_gene_to_data, flank_size, min_seq_len):
    """Cut a `flank_size`-padded window around each row's sense start; format as RNA.

    Rows with a missing gene, or a flanked window shorter than `min_seq_len`, are
    dropped. Returns (rna_seqs, sense_len_map, sense_start_map), all keyed by the
    positional `_temp_id`: the (id, RNA sequence) pairs, the sense length, and the
    sense start relative to the flanked window.
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

    For each rna_id, slide an `access_size` window over the raccess energies
    (window_access_energies) and average the per-seed columns into a single
    `avg_access`. Seeds larger than `access_size`, and sequences shorter than it,
    contribute nothing. Returns a DataFrame ["rna_id", "avg_access"] in per-rna_id
    position order.
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
        df = window_access_energies(access_res, rna_size, access_size, eligible_seeds)
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


def _run_batches(process_batch, valid_df, batch_size, n_jobs):
    """Slice valid_df into batches and run `process_batch` over them.

    Threaded when n_jobs > 1 — raccess runs as a subprocess, so the GIL doesn't
    serialize the workers.
    """
    batches = [valid_df.iloc[i : i + batch_size] for i in range(0, len(valid_df), batch_size)]
    if n_jobs > 1:
        with ThreadPoolExecutor(max_workers=min(n_jobs, len(batches))) as pool:
            return list(pool.map(process_batch, batches))
    return [process_batch(b) for b in batches]


def populate_sense_accessibility_multi(
    aso_dataframe,
    gene_to_data,
    configs,
    batch_size=1000,
    access_win_size=ACCESS_WIN_SIZE,
    n_jobs=1,
):
    """Populate one accessibility feature column per config.

    Configs sharing the same ``flank`` are processed in a single pass: the
    flank-padded window is cut once per batch and ONE raccess call covers the
    union of seed sizes; each config then takes its own access-window slice.
    Configs with different flanks are processed in separate per-flank passes
    (one raccess call per (batch, flank)). The returned ``feature_names`` list
    mirrors the order of ``configs``.
    """
    if not configs:
        return aso_dataframe.copy(), []

    # Group by flank, preserving the input order so feature_names mirrors it.
    flank_groups = defaultdict(list)
    for idx, cfg in enumerate(configs):
        flank_groups[cfg["flank"]].append((idx, cfg))

    df_out = aso_dataframe.copy()
    name_by_idx = {}
    for indexed_cfgs in flank_groups.values():
        orig_idxs, group_cfgs = zip(*indexed_cfgs)
        group_names = _populate_single_flank_accessibility(
            df_out,
            gene_to_data,
            list(group_cfgs),
            batch_size=batch_size,
            access_win_size=access_win_size,
            n_jobs=n_jobs,
        )
        for idx, fn in zip(orig_idxs, group_names):
            name_by_idx[idx] = fn

    feature_names = [name_by_idx[i] for i in range(len(configs))]
    return df_out, feature_names


def _populate_single_flank_accessibility(
    df_out,
    gene_to_data,
    configs,
    batch_size,
    access_win_size,
    n_jobs,
):
    """Add per-config accessibility columns to ``df_out`` in place; all configs must share ``flank``."""
    flank_size = configs[0]["flank"]
    union_seeds = sorted({s for c in configs for s in c["seeds"]})
    min_access = min(c["access"] for c in configs)
    feature_names = [_build_feature_name(c["flank"], c["access"], c["seeds"]) for c in configs]

    df_out["_temp_id"] = range(len(df_out))

    valid_mask = df_out[SENSE_START] != -1
    if not valid_mask.any():
        for fn in feature_names:
            df_out[fn] = np.nan
        df_out.drop(columns=["_temp_id"], inplace=True)
        return feature_names

    valid_df = df_out.loc[valid_mask, ["_temp_id", SENSE_START, SENSE_LENGTH, CANONICAL_GENE]].copy()
    lightweight_gene_to_data = _lightweight_gene_to_data(valid_df[CANONICAL_GENE].unique(), gene_to_data)
    raccess_exe = find_raccess()

    def _process_batch(batch_rows):
        rna_seqs, sense_len_map, sense_start_map = _build_batch_rna_seqs(
            batch_rows, lightweight_gene_to_data, flank_size, min_seq_len=min_access
        )
        if not rna_seqs:
            return {fn: pd.DataFrame(columns=["_temp_id", fn]) for fn in feature_names}

        # ONE raccess call for this batch over the union of all configs' seeds.
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

    batch_results = _run_batches(_process_batch, valid_df, batch_size, n_jobs)
    for fn in feature_names:
        _assign_feature_column(df_out, [br[fn] for br in batch_results], fn)
    df_out.drop(columns=["_temp_id"], inplace=True)
    return feature_names


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
    """Populate ONE accessibility feature column (single-config wrapper over _multi)."""
    if seed_sizes is None:
        seed_sizes = SEED_SIZES
    seeds_list = list(seed_sizes) if isinstance(seed_sizes, (list, tuple)) else [seed_sizes]
    config = {"flank": flank_size, "access": access_size, "seeds": seeds_list}

    df_out, feature_names = populate_sense_accessibility_multi(
        aso_dataframe,
        gene_to_data,
        [config],
        batch_size=batch_size,
        access_win_size=access_win_size,
        n_jobs=n_jobs,
    )
    return df_out, feature_names[0]
