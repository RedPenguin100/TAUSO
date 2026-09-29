"""Assemble everything Figure 2's panels read, into data/.

The panels never touch the model, the feature store or the training data: they read the three
files written here. That keeps a panel reproducible from the folder alone, and makes it explicit
when the figure is redrawn from a different model -- the inputs change in one visible commit
rather than silently under the plotting code.

    python notebooks/figures/figure2/build_data.py [--k 5] [--min-screen 25]

Writes:
    data/scores.parquet          one row per ASO in the comparison subset: labels, TAUSO's
                                 prediction and every competitor score, already sign-oriented
    data/screening_effort.csv    per screen, the random-screening fold each method saves
    data/parsimony.csv           the top-K feature sweep (produced separately; see README)
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))
from notebooks.models import common  # noqa: E402

from tauso.data.consts import CANONICAL_GENE_NAME, CELL_LINE, CHEMICAL_PATTERN, INHIBITION_PERCENT  # noqa: E402
from tauso.data.data import get_data_dir  # noqa: E402
from tauso.inference.scoring import load_model  # noqa: E402

DATA = HERE / "data"
IDX, CID = "index_oligo", "custom_id"
STRICT = ("MMMMMddddddddddMMMMM", "CCCddddddddddCCC")

# label -> the shard column holding that tool's score
COMPETITORS = {
    "OligoAI": "oligo_ai_score",
    "OligoWalk": "OW_Overall",
    "OligoWalk·Tm": "OW_Tm",
    "OligoWalk·intra": "OW_Intra_Oligo",
    "OligoWalk·duplex": "OW_Duplex",
    "miRanda": "miranda_score",
    "miRanda·E": "miranda_energy",
    "sfold": "sfold_accessibility",
    "PFRED": "PFRED_PLS",
    "PFRED·SVM": "PFRED_SVM",
}


def _competitor_shards():
    """Every competitor score keyed on index_oligo, from the feature store's competition shards."""
    folder = Path(get_data_dir()) / "features" / "oligo" / "_competitors_v15"
    if not folder.is_dir():
        sys.exit(f"{folder} not found. Fetch it with `tauso setup-features --include-competition`.")
    merged = None
    for label, column in COMPETITORS.items():
        parquet, csv = folder / f"{column}.parquet", folder / f"{column}.csv"
        if parquet.exists():
            shard = pd.read_parquet(parquet)
        elif csv.exists():
            shard = pd.read_csv(csv)
        else:
            print(f"  ! {label}: no shard, skipping")
            continue
        shard = shard[[IDX, column]].rename(columns={column: label})
        merged = shard if merged is None else merged.merge(shard, on=IDX, how="outer")
    return merged


def _median_spearman(frame, column, group):
    """Median Spearman of a score against measured inhibition, over screens of at least 3 ASOs."""
    from scipy.stats import spearmanr

    sub = frame[[column, INHIBITION_PERCENT, group]].dropna()
    values = [
        spearmanr(part[column], part[INHIBITION_PERCENT]).correlation
        for _, part in sub.groupby(group)
        if part[column].nunique() > 1 and len(part) >= 3
    ]
    values = [v for v in values if np.isfinite(v)]
    return float(np.median(values)) if values else np.nan


def _orient(frame):
    """Flip a score stored the other way round, so that larger always argues for more knockdown.

    A tool is flipped only when it is anti-correlated in both groupings, which distinguishes a
    reversed convention from a tool that simply ranks badly in one of them.
    """
    signs = {}
    for label in [c for c in COMPETITORS if c in frame.columns]:
        per_screen = _median_spearman(frame, label, CID)
        per_cohort = _median_spearman(frame, label, "cohort")
        backwards = np.isfinite(per_screen) and np.isfinite(per_cohort) and per_screen < 0 and per_cohort < 0
        signs[label] = -1.0 if backwards else 1.0
    return signs


def _expected_best_of(knockdowns, draws, k, boot, rng):
    """Mean knockdown of the best k among `draws` ASOs picked at random, averaged over `boot` draws."""
    picked = np.argsort(rng.random((boot, len(knockdowns))), axis=1)[:, :draws]
    return np.sort(knockdowns[picked], axis=1)[:, -k:].mean()


def _fold_saved(knockdowns, order, k, boot, rng):
    """How many ASOs a random screen must test to match this method's top k, divided by k."""
    target = knockdowns[order[:k]].mean()
    low, high = k, len(knockdowns)
    while low < high:
        middle = (low + high) // 2
        if _expected_best_of(knockdowns, middle, k, boot, rng) >= target:
            high = middle
        else:
            low = middle + 1
    return low / k


def screening_effort(frame, k, min_screen, boot, rng):
    """Per screen, the factor by which each method shrinks a random screen matching its top k."""
    methods = {"TAUSO": "TAUSO", "OligoAI": "OligoAI", "random": None}
    folds = {name: [] for name in methods}
    for _, screen in frame.groupby(CID):
        knockdowns = screen[INHIBITION_PERCENT].to_numpy(np.float64)
        if len(screen) < min_screen or np.unique(knockdowns).size < 2:
            continue
        for name, column in methods.items():
            if column is None:
                order = rng.permutation(len(knockdowns))
            else:
                order = np.argsort(-screen[column].to_numpy(np.float64))
            folds[name].append(_fold_saved(knockdowns, order, k, boot, rng))
    return pd.DataFrame(folds)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--k", type=int, default=5, help="shortlist size for the screening-effort metric")
    ap.add_argument("--min-screen", type=int, default=25, help="smallest screen the metric uses")
    ap.add_argument("--boot", type=int, default=1500, help="random draws per screen size")
    args = ap.parse_args()

    booster, features = load_model()
    frame, _ = common.load_dataset()
    _, test = common.split(frame)
    print(f"held-out test: {len(test):,} ASOs, {len(features)} features")

    scored = test[[IDX, INHIBITION_PERCENT, CID, CANONICAL_GENE_NAME, CELL_LINE, CHEMICAL_PATTERN]].copy()
    scored["TAUSO"] = common.predict(booster, test, features)
    scored = scored.merge(_competitor_shards(), on=IDX, how="left")

    # The shared subset: canonical gapmer architectures, and only rows every tool scored, so no
    # method is credited or penalised for the rows it happens to cover.
    scored = scored[scored[CHEMICAL_PATTERN].isin(STRICT) & scored["OligoAI"].notna()].reset_index(drop=True)
    scored["cohort"] = scored[CANONICAL_GENE_NAME].astype(str) + "|" + scored[CELL_LINE].astype(str)

    signs = _orient(scored)
    for label, sign in signs.items():
        scored[label] = sign * scored[label]
    flipped = [label for label, sign in signs.items() if sign < 0]
    print(f"subset: {len(scored):,} ASOs, {scored[CID].nunique()} screens, {scored['cohort'].nunique()} cohorts")
    print(f"sign-flipped (stored the other way round): {flipped or 'none'}")

    DATA.mkdir(parents=True, exist_ok=True)
    scored.to_parquet(DATA / "scores.parquet", index=False)
    print(f"wrote {DATA / 'scores.parquet'}")

    effort = screening_effort(scored, args.k, args.min_screen, args.boot, np.random.default_rng(0))
    effort.to_csv(DATA / "screening_effort.csv", index=False)
    geomean = {c: float(np.exp(np.log(effort[c]).mean())) for c in effort}
    print(f"wrote {DATA / 'screening_effort.csv'} ({len(effort)} screens, k={args.k})")
    print("  " + "  ".join(f"{name} {value:.2f}x" for name, value in geomean.items()))


if __name__ == "__main__":
    main()
