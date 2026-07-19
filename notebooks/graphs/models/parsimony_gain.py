"""Parsimony sweep: held-out within-experiment Spearman vs number of top-K features.

Features are ranked by GAIN importance taken from the shipped booster. For each K the champion
config is retrained on train+val using only the top-K features and scored once on the held-out
test, restricted to the same comparison subset Fig 2 uses (strict gapmers, rows OligoAI also
scored) so the competitor reference lines are directly comparable.

Writes parsimony_gain_test.csv incrementally (one row per K) so a long run is never lost.

    python parsimony_gain.py [--ks 485,450,...] [--seed 1]
"""
import argparse
import json
import sys
import time
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from notebooks.models import common                 # noqa: E402
from tauso.inference import load_model              # noqa: E402
import _modeldata as M                              # noqa: E402

OUT = Path(__file__).resolve().parent / "parsimony_gain_test.csv"
BEST = json.loads((REPO / "notebooks/models/best_models_parameters.json").read_text())["clean_exp_exp_med"]
KS = [485, 450, 400, 350, 300, 275, 250, 225, 200, 180, 160, 150, 140, 120,
      100, 90, 80, 70, 60, 50, 40, 35, 30, 25, 20]


def gain_ranked_features(model, features):
    """Model features ordered by descending gain; split-less features (gain 0) keep a stable tail."""
    gain = model.get_score(importance_type="gain")
    return sorted(features, key=lambda f: -gain.get(f, 0.0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ks", default=None, help="comma-separated K values (default: the selector grid)")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    ks = [int(k) for k in args.ks.split(",")] if args.ks else KS

    df, _ = common.load_dataset()
    trv, test = common.split(df)
    model, features = load_model()
    ranked = gain_ranked_features(model, features)
    print(f"train+val={len(trv)}  test={len(test)}  features={len(features)}", flush=True)

    sub, present = M.load()                      # the Fig 2 comparison subset
    sign = M.orient(sub, present)
    refs = {lbl: M.med(sub, col, M.CID) * sign[lbl] for lbl, col in present.items() if lbl != "TAUSO"}
    print(f"eval subset n={len(sub)}  refs: " +
          " ".join(f"{k}={v:.3f}" for k, v in sorted(refs.items(), key=lambda kv: -kv[1])), flush=True)

    params, rounds, variant = BEST["params"], BEST["num_boost_round"], BEST["variant"]
    rows = []
    for i, k in enumerate(ks, 1):
        t0 = time.time()
        feats = ranked[:k]
        m = common.train(trv, feats, variant, params, rounds, seed=args.seed)
        pred = common.predict(m, test, feats)
        scored = sub.merge(pd.DataFrame({M.IDX: test["index_oligo"].to_numpy(), "score": pred}),
                           on=M.IDX, how="left")
        em = M.med(scored, "score", M.CID)
        gx = M.med(scored, "score", "gxc")
        rows.append({"K": k, "exp_med": em, "gxc_med": gx, "seconds": round(time.time() - t0, 1)})
        pd.DataFrame(rows).to_csv(OUT, index=False)          # incremental save
        print(f"[{i}/{len(ks)}] K={k:4d}  exp_med={em:.4f}  gxc_med={gx:.4f}  ({time.time()-t0:.0f}s)", flush=True)

    pd.DataFrame([{"label": k, "exp_med": v} for k, v in refs.items()]).to_csv(
        OUT.with_name("parsimony_gain_refs.csv"), index=False)
    print(f"\nwrote {OUT} and {OUT.with_name('parsimony_gain_refs.csv')}", flush=True)


if __name__ == "__main__":
    main()
