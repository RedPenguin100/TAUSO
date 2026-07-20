"""Screening-effort reduction, TAUSO vs OligoAI, on the shared held-out test.

OligoAI's metric (KCNT2 validation): to match the knockdown quality of a tool's top-k picks, how many
ASOs must a random screen test (keeping its best k)? fold = N_random / k. We compute the SAME mechanic
per test screen for both tools on the identical OligoAI split, then aggregate.

    python screening_effort.py [--k 5] [--min 25] [--boot 1500]
"""
import argparse
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(REPO / "src")); sys.path.insert(0, str(Path(__file__).resolve().parent))
import _modeldata as M   # noqa: E402


def E_bestk(kd, N, k, B, rng):
    """E[ mean of the best k among N random draws without replacement ] over B bootstraps."""
    perm = np.argsort(rng.random((B, len(kd))), axis=1)[:, :N]
    samp = np.sort(kd[perm], axis=1)[:, -k:]
    return samp.mean()


def fold_for(kd, order, k, B, rng):
    """Smallest N with E[best-k-of-N random] >= mean KD of the top-k *predicted*, / k."""
    target = kd[order[:k]].mean()
    n = len(kd)
    lo, hi = k, n
    while lo < hi:
        mid = (lo + hi) // 2
        if E_bestk(kd, mid, k, B, rng) >= target:
            hi = mid
        else:
            lo = mid + 1
    return lo / k


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--min", type=int, default=25)
    ap.add_argument("--boot", type=int, default=1500)
    args = ap.parse_args()
    rng = np.random.default_rng(0)

    d, present = M.load()
    sign = M.orient(d, present)
    sizes = d.groupby(M.CID).size()
    print(f"screens: {d[M.CID].nunique()}  | with >= {args.min} ASOs: {(sizes >= args.min).sum()}  "
          f"| size median {int(sizes.median())} max {int(sizes.max())}", flush=True)

    tools = {"TAUSO": ("tauso", 1.0), "OligoAI": ("oligo_ai_score", sign["OligoAI"]),
             "random": (None, 1.0)}
    folds = {t: [] for t in tools}
    for cid, s in d.groupby(M.CID):
        if len(s) < args.min:
            continue
        kd = s[M.INH].to_numpy(np.float64)
        if np.unique(kd).size < 2:
            continue
        for t, (col, sg) in tools.items():
            order = rng.permutation(len(kd)) if col is None else np.argsort(-(sg * s[col].to_numpy(np.float64)))
            folds[t].append(fold_for(kd, order, args.k, args.boot, rng))

    print(f"\nscreening-effort reduction (k={args.k}, {len(folds['TAUSO'])} screens):")
    print(f"{'method':>8} {'geomean':>9} {'median':>8} {'IQR':>15}")
    for t in tools:
        f = np.array(folds[t])
        gm = np.exp(np.mean(np.log(f)))
        q1, md, q3 = np.percentile(f, [25, 50, 75])
        print(f"{t:>8} {gm:>8.2f}x {md:>7.2f}x  {q1:.2f}-{q3:.2f}", flush=True)
    import pandas as pd
    pd.DataFrame(folds).to_csv(Path(__file__).resolve().parent / "screening_effort.csv", index=False)


if __name__ == "__main__":
    main()
