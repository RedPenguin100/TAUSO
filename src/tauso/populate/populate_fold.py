import logging
from collections import defaultdict

import numpy as np
import pandas as pd

from ..data.consts import CANONICAL_GENE_NAME, STRUCTURE_SENSE_LENGTH, STRUCTURE_SENSE_START
from ..features.fold.vienna_access import calculate_avg_access_per_setting
from ..features.fold.vienna_fold import calculate_avg_mfe_per_setting, calculate_end_mfe
from ..parallel_utils import make_apply_fn

logger = logging.getLogger(__name__)


def validate_cols_in_df(df, cols):
    missing_cols = [c for c in cols if c not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing columns in DataFrame: {missing_cols}")


# MFE grid: (flank, window, step). Six configs spanning three flank scales with
# two window sizes each, covering tight (sub-target), local (stem-loop), and
# mid-regional structural contexts. Window size is the dominant axis; the paired
# windows at each flank explore the "tighter vs slightly wider" trade-off
# (e.g. catching a small hairpin only inside the wider window).
DEFAULT_SETTINGS = [
    (30, 25, 4),  # very tight: ASO-target-sized window
    (30, 40, 5),  # tight: sub-target stems with small flanking context
    (60, 55, 5),  # local: short stems inside a mid-flank context
    (60, 70, 5),  # local: canonical stem-loop scale
    (120, 100, 7),  # mid-regional
    (120, 150, 10),  # wider-regional
]


# Target sites this close together on one gene share a single fold. A wider chunk serves
# more sites per fold but costs more to fold; around 500nt the two balance out.
MFE_CHUNK_SIZE = 500

# The two ASO ends contact different target sub-windows and are not interchangeable:
# hybridisation nucleates at whichever terminus makes first contact, while RNase H1 then
# cleaves in a fixed register from the ASO-3' end. Both readouts come off folds the grids
# above already do, so the four end features cost no additional folding.
END_LEN = 6
MFE_END_SETTING = (30, 40, 5)

FOLD_REGION_START = "_mfe_fold_start"
FOLD_REGION_END = "_mfe_fold_end"


def _lightweight_gene_to_data(genes, gene_to_data):
    """gene -> mRNA string, so worker threads don't pickle the heavy gene_to_data."""
    return {gene: str(gene_to_data[gene].full_mrna) for gene in genes if gene in gene_to_data}


def _fold_work_frame(df, lightweight_gene_to_data, widest_flank, chunk_size):
    """Just the columns the MFE apply reads, plus the region of the gene to fold per row.

    Target sites on the same gene whose flanked cuts all fit within `chunk_size` are given
    the same region, so one fold serves all of them. Rows come back in region order, so
    that fold is still cached when the next row asks for it.
    """
    work = df[[CANONICAL_GENE_NAME, STRUCTURE_SENSE_START, STRUCTURE_SENSE_LENGTH]].copy()

    sites_by_gene = defaultdict(list)
    for label, gene, start, length in work.itertuples(name=None):
        if gene in lightweight_gene_to_data and start != -1:
            sites_by_gene[gene].append((start, length, label))

    regions = {}
    for gene, sites in sites_by_gene.items():
        mrna_len = len(lightweight_gene_to_data[gene])
        sites.sort()
        i = 0
        while i < len(sites):
            region_start = max(0, sites[i][0] - widest_flank)
            region_end = 0
            batch = []
            while i < len(sites):
                start, length, label = sites[i]
                site_end = min(mrna_len, start + length + widest_flank)
                if batch and site_end - region_start > chunk_size:
                    break
                region_end = max(region_end, site_end)
                batch.append(label)
                i += 1
            for label in batch:
                regions[label] = (region_start, region_end)

    work[FOLD_REGION_START] = [regions.get(label, (0, 0))[0] for label in work.index]
    work[FOLD_REGION_END] = [regions.get(label, (0, 0))[1] for label in work.index]
    return work.sort_values([CANONICAL_GENE_NAME, FOLD_REGION_START], kind="stable")


def mfe_feature_name(flank, window_size, step):
    """Stable column name for one MFE setting."""
    return f"fold_mfe_win{window_size}_flank{flank}_step{step}"


def populate_mfe_features(df, gene_to_data, n_jobs=1, verbose=False, settings=None):
    if settings is None:
        settings = DEFAULT_SETTINGS

    required_cols = [CANONICAL_GENE_NAME, STRUCTURE_SENSE_START, STRUCTURE_SENSE_LENGTH]
    validate_cols_in_df(df, required_cols)

    lightweight_gene_to_data = _lightweight_gene_to_data(df[CANONICAL_GENE_NAME].dropna().unique(), gene_to_data)

    # All settings are handled in a single apply. Every setting's sub-sequence cut sits
    # inside the widest one, so one pass per row both shares the folding across settings
    # and pays the parallel-dispatch overhead exactly once.
    feature_names = [mfe_feature_name(f, w, s) for f, w, s in settings]
    end_names = [f"fold_mfe_{end}" for end in ("aso5end", "aso3end")] if MFE_END_SETTING in settings else []

    widest_flank = max(flank for flank, _, _ in settings)
    work = _fold_work_frame(df, lightweight_gene_to_data, widest_flank, MFE_CHUNK_SIZE)

    def _process_row(row):
        gene_name = row[CANONICAL_GENE_NAME]
        global_start = row[STRUCTURE_SENSE_START]
        sense_len = row[STRUCTURE_SENSE_LENGTH]

        out = {name: np.nan for name in feature_names + end_names}
        if gene_name not in lightweight_gene_to_data or global_start == -1:
            return out
        full_mrna = lightweight_gene_to_data[gene_name]

        per_setting = calculate_avg_mfe_per_setting(
            full_mrna,
            global_start,
            sense_len,
            settings,
            fold_region=(row[FOLD_REGION_START], row[FOLD_REGION_END]),
        )
        for (flank_size, window_size, step), value in per_setting.items():
            out[mfe_feature_name(flank_size, window_size, step)] = value
        if end_names:
            ends = calculate_end_mfe(
                full_mrna,
                global_start,
                sense_len,
                *MFE_END_SETTING,
                END_LEN,
                fold_region=(row[FOLD_REGION_START], row[FOLD_REGION_END]),
            )
            for end, value in ends.items():
                out[f"fold_mfe_{end}"] = value
        return out

    apply_fn = make_apply_fn(work, n_jobs=n_jobs, progress_bar=verbose, verbose=2 if verbose else 0)
    results = apply_fn(_process_row, axis=1)
    results_df = pd.DataFrame(list(results), index=work.index)
    feature_names = feature_names + end_names
    for name in feature_names:
        df[name] = results_df[name]
    return df, feature_names


# Accessibility grid: (flank, max_bp_span, open_len, anchor). Every setting shares one
# fold of the flank-60 cut, so opening lengths and anchorings are free once it is done
# and the grid takes the whole usable range of both. Flank is the axis that costs, and
# 60 is where the curve flattens. A span of None leaves the 140nt cut free to pair across
# itself -- a number wider than the cut would read as unconstrained while silently biting
# once targets outgrew it.
DEFAULT_ACCESS_SETTINGS = [
    # `a5` and `a3` sweep the whole target, one window per position.
    (60, None, open_len, anchor)
    for anchor in ("a5", "a3")
    for open_len in (4, 6, 8, 10, 13, 16, 20, 26, 32)
] + [
    # `aso5end` and `aso3end` take the single window flush with each end of it.
    (60, None, END_LEN, anchor)
    for anchor in ("aso5end", "aso3end")
]


def access_feature_name(flank, max_bp_span, open_len, anchor):
    """Stable column name for one accessibility setting; None span is written `sinf`."""
    span = "inf" if max_bp_span is None else max_bp_span
    return f"access_f{flank}_s{span}_u{open_len}_{anchor}"


def populate_access_features(df, gene_to_data, n_jobs=1, verbose=False, settings=None):
    """Add one accessibility column per (flank, max_bp_span, open_len) setting.

    Rows with an unknown gene, or a sense start of -1, come out NaN.
    """
    if settings is None:
        settings = DEFAULT_ACCESS_SETTINGS

    required_cols = [CANONICAL_GENE_NAME, STRUCTURE_SENSE_START, STRUCTURE_SENSE_LENGTH]
    validate_cols_in_df(df, required_cols)

    lightweight_gene_to_data = _lightweight_gene_to_data(df[CANONICAL_GENE_NAME].dropna().unique(), gene_to_data)

    feature_names = [access_feature_name(*setting) for setting in settings]

    def _process_row(row):
        gene_name = row[CANONICAL_GENE_NAME]
        sense_start = row[STRUCTURE_SENSE_START]
        sense_len = row[STRUCTURE_SENSE_LENGTH]

        out = {name: np.nan for name in feature_names}
        # If we can't identify the target RNA
        if gene_name not in lightweight_gene_to_data or sense_start == -1:
            return out
        full_mrna = lightweight_gene_to_data[gene_name]

        per_setting = calculate_avg_access_per_setting(full_mrna, sense_start, sense_len, settings)
        for setting, value in per_setting.items():
            out[access_feature_name(*setting)] = value
        return out

    # ViennaRNA holds the GIL for the duration of a fold, so this has to be
    # process parallelism -- threads measured identical at n_jobs 1, 8 and 32.
    apply_fn = make_apply_fn(df, n_jobs=n_jobs, progress_bar=verbose, verbose=2 if verbose else 0)
    results = apply_fn(_process_row, axis=1)
    results_df = pd.DataFrame(list(results), index=df.index)
    for name in feature_names:
        df[name] = results_df[name]
    return df, feature_names
