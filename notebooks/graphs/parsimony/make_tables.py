"""Builds the per-rung result table and the drop-order table from curve.csv, test_curve.csv and feats/."""
import json
import re
from pathlib import Path

import pandas as pd

L = Path(__file__).resolve().parent

cv = pd.read_csv(L / "curve.csv")
te = pd.read_csv(L / "test_curve.csv")
cols = ["rmse", "mae", "exp_med", "exp_mean", "gxc_med", "gxc_mean"]
tcols = ["rmse", "mae", "exp_med", "exp_mean", "gxc_med", "gxc_mean", "top5", "p5", "p10", "globP"]
table = cv[["n"] + [f"cv_{c}" for c in cols]].merge(
    te[["n"] + [f"test_{c}" for c in tcols]], on="n")
table.sort_values("n", ascending=False).round(4).to_csv(L / "results_by_n.csv", index=False)

sets = {int(re.search(r"feats_(\d+)", p.name).group(1)): json.loads(p.read_text())
        for p in (L / "feats").glob("feats_*.json")}
ns = sorted(sets, reverse=True)
rows = []
for big, small in zip(ns, ns[1:]):
    gone = sorted(set(sets[big]) - set(sets[small]))
    rows += [{"feature": f, "last_n_present": big, "first_n_absent": small,
              "dropped_together_with": len(gone) - 1} for f in gone]
rows += [{"feature": f, "last_n_present": ns[-1], "first_n_absent": None,
          "dropped_together_with": None} for f in sorted(sets[ns[-1]])]
drop = pd.DataFrame(rows).sort_values(["last_n_present", "feature"], ascending=[False, True])
drop[["first_n_absent", "dropped_together_with"]] = drop[["first_n_absent", "dropped_together_with"]].astype("Int64")
drop.to_csv(L / "drop_order.csv", index=False)
print(f"{len(table)} rungs, {len(drop)} features ({drop.first_n_absent.notna().sum()} dropped, "
      f"{drop.first_n_absent.isna().sum()} kept at n={ns[-1]})")
