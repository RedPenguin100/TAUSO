import logging
from collections import defaultdict

import numpy as np
import pandas as pd

from ..data.consts import CANONICAL_GENE_NAME, STRUCTURE_SENSE_LENGTH, STRUCTURE_SENSE_START
from ..features.fold.vienna_access import calculate_avg_access_per_setting
from ..features.fold.vienna_fold import calculate_avg_mfe_per_setting, calculate_end_mfe
from ..genome.read_human_genome import get_gene_to_data_subset
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

# How many nucleotides accessibility holds unpaired for the `aso5end`/`aso3end` anchors and
# the two `std` features in the grid below.
END_LEN = 6

# How many terminal target positions each MFE end feature averages over. The two ASO ends
# contact different target sub-windows and are not interchangeable: hybridisation nucleates
# at whichever terminus makes first contact, while RNase H1 then cleaves in a fixed register
# from the ASO-3' end.
TERMINAL_MFE_DEFAULT = 6
MFE_END_SETTING = (30, 40, 5)

FOLD_REGION_START = "_mfe_fold_start"
FOLD_REGION_END = "_mfe_fold_end"


def _plan_fold_regions(df, gene_to_premrna, widest_flank, chunk_size):
    """Decide which region of the gene each row folds, and order the rows to match.

    Target sites on the same gene whose flanked cuts all fit within `chunk_size` are given
    the same region, so one fold serves all of them. The returned rows are sorted by that
    region, which is what makes the sharing pay: `_folded` caches only 8 folds, so out of
    order the fold would be evicted before the next row in its region asked for it.

    Carries only the three columns the apply reads, so workers do not pickle the rest.
    """
    planned = df[[CANONICAL_GENE_NAME, STRUCTURE_SENSE_START, STRUCTURE_SENSE_LENGTH]].copy()

    sites_by_gene = defaultdict(list)
    for label, gene, start, length in planned.itertuples(name=None):
        if gene in gene_to_premrna and start != -1:
            sites_by_gene[gene].append((start, length, label))

    regions = {}
    for gene, sites in sites_by_gene.items():
        mrna_len = len(gene_to_premrna[gene])
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

    planned[FOLD_REGION_START] = [regions.get(label, (0, 0))[0] for label in planned.index]
    planned[FOLD_REGION_END] = [regions.get(label, (0, 0))[1] for label in planned.index]
    return planned.sort_values([CANONICAL_GENE_NAME, FOLD_REGION_START], kind="stable")


def _populate_per_row(df, gene_to_mrna, feature_names, compute, n_jobs, verbose, planned=None):
    """Fill `feature_names` from `compute(full_mrna, sense_start, sense_len, row)` per row.

    Rows whose gene is unknown, or whose sense start is -1, come out NaN without calling
    `compute`. `planned` is the frame the apply runs over when it differs from `df`: the MFE
    path hands in one sorted by fold region so neighbouring sites share a fold, while the
    accessibility path folds each row's own cut and needs no ordering.

    ViennaRNA holds the GIL for the duration of a fold, so this has to be process
    parallelism -- threads measured identical at n_jobs 1, 8 and 32.
    """
    frame = df if planned is None else planned

    def _process_row(row):
        out = {name: np.nan for name in feature_names}
        gene_name = row[CANONICAL_GENE_NAME]
        sense_start = row[STRUCTURE_SENSE_START]
        if gene_name not in gene_to_mrna or sense_start == -1:
            return out
        out.update(compute(gene_to_mrna[gene_name], sense_start, row[STRUCTURE_SENSE_LENGTH], row))
        return out

    apply_fn = make_apply_fn(frame, n_jobs=n_jobs, progress_bar=verbose, verbose=2 if verbose else 0)
    results_df = pd.DataFrame(list(apply_fn(_process_row, axis=1)), index=frame.index)
    for name in feature_names:
        df[name] = results_df[name]
    return df, feature_names


def mfe_feature_name(flank, window_size, step):
    """Stable column name for one MFE setting."""
    return f"fold_mfe_win{window_size}_flank{flank}_step{step}"


def populate_mfe_features(df, gene_to_data, n_jobs=1, verbose=False, settings=None):
    if settings is None:
        settings = DEFAULT_SETTINGS

    required_cols = [CANONICAL_GENE_NAME, STRUCTURE_SENSE_START, STRUCTURE_SENSE_LENGTH]
    validate_cols_in_df(df, required_cols)

    gene_to_mrna = get_gene_to_data_subset(df[CANONICAL_GENE_NAME].dropna().unique(), gene_to_data)

    # All settings are handled in a single apply. Every setting's sub-sequence cut sits
    # inside the widest one, so one pass per row both shares the folding across settings
    # and pays the parallel-dispatch overhead exactly once.
    end_names = [f"fold_mfe_{k}" for k in ("aso5end", "aso3end", "std")] if MFE_END_SETTING in settings else []
    feature_names = [mfe_feature_name(f, w, s) for f, w, s in settings] + end_names

    widest_flank = max(flank for flank, _, _ in settings)
    planned = _plan_fold_regions(df, gene_to_mrna, widest_flank, MFE_CHUNK_SIZE)

    def compute(full_mrna, sense_start, sense_len, row):
        fold_region = (row[FOLD_REGION_START], row[FOLD_REGION_END])
        out = {
            mfe_feature_name(*setting): value
            for setting, value in calculate_avg_mfe_per_setting(
                full_mrna, sense_start, sense_len, settings, fold_region
            ).items()
        }
        if end_names:
            ends = calculate_end_mfe(
                full_mrna, sense_start, sense_len, *MFE_END_SETTING, TERMINAL_MFE_DEFAULT, fold_region
            )
            out.update({f"fold_mfe_{end}": value for end, value in ends.items()})
        return out

    return _populate_per_row(df, gene_to_mrna, feature_names, compute, n_jobs, verbose, planned=planned)


# Accessibility grid: (flank, max_bp_span, open_len, anchor, reducer). Every setting
# shares one fold of the flank-60 cut, so opening lengths, anchorings and reducers are
# free once it is done and the grid takes the whole usable range of them. Flank is the
# axis that costs, and 60 is where the curve flattens. A span of None leaves the 140nt
# cut free to pair across itself -- a number wider than the cut would read as
# unconstrained while silently biting once targets outgrew it.
ACCESS_OPEN_LENS = (4, 6, 8, 10, 13, 16, 20, 26, 32)

DEFAULT_ACCESS_SETTINGS = (
    [
        # `a5` and `a3` sweep the whole target, one window per position.
        (60, None, open_len, anchor, "mean")
        for anchor in ("a5", "a3")
        for open_len in ACCESS_OPEN_LENS
    ]
    + [
        # `aso5end` and `aso3end` take the single window flush with each end of it. One
        # window has no spread, so these are mean only.
        (60, None, END_LEN, anchor, "mean")
        for anchor in ("aso5end", "aso3end")
    ]
    + [
        # How evenly open the target is, rather than how open. Taken at one opening length
        # under both anchorings: spreads at neighbouring opening lengths measure nearly the
        # same thing, and carrying all nine scored worse than carrying one.
        (60, None, END_LEN, anchor, "std")
        for anchor in ("a5", "a3")
    ]
)


def access_feature_name(flank, max_bp_span, open_len, anchor, reducer):
    """Stable column name for one accessibility setting; None span is written `sinf`.

    The mean is left unmarked in the name so that adding a reducer to the grid does not
    rename columns that already exist.
    """
    span = "inf" if max_bp_span is None else max_bp_span
    name = f"access_f{flank}_s{span}_u{open_len}_{anchor}"
    if reducer == "mean":
        return name
    if reducer == "std":
        return f"{name}_std"
    raise ValueError(f"unknown reducer {reducer!r}; expected one of mean, std")


def populate_access_features(df, gene_to_data, n_jobs=1, verbose=False, settings=None):
    """Add one accessibility column per (flank, max_bp_span, open_len) setting."""
    if settings is None:
        settings = DEFAULT_ACCESS_SETTINGS

    required_cols = [CANONICAL_GENE_NAME, STRUCTURE_SENSE_START, STRUCTURE_SENSE_LENGTH]
    validate_cols_in_df(df, required_cols)

    gene_to_mrna = get_gene_to_data_subset(df[CANONICAL_GENE_NAME].dropna().unique(), gene_to_data)

    feature_names = [access_feature_name(*setting) for setting in settings]

    def compute(full_mrna, sense_start, sense_len, row):
        per_setting = calculate_avg_access_per_setting(full_mrna, sense_start, sense_len, settings)
        return {access_feature_name(*setting): value for setting, value in per_setting.items()}

    return _populate_per_row(df, gene_to_mrna, feature_names, compute, n_jobs, verbose)
