"""The competitor reference rows: every published tool's score on the same held-out test rows.

The competitor scores are not computed here. They ship as shards on the feature record and are
fetched and md5-checked by `ensure_competition_cache`, the same way the rest of the repo gets
them; this reads those columns and evaluates each one exactly like every other row in the
comparison. `COMPETITION` (notebooks/features/feature_extraction) is the list of them.

Sanity floor: OligoAI is worth about 0.40 exp_med on this data. A value near zero means the
`index_oligo` join misaligned, not that the tool is bad -- see the compare-to-oligoai notes.

Run (needs the tauso package and TAUSO_DATA_DIR):
  python notebooks/competitors/oligogym/score_competitors.py
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
sys.path[:0] = [str(REPO_ROOT), str(REPO_ROOT / "src")]

from notebooks.features.feature_extraction import COMPETITION  # noqa: E402
from notebooks.models import common  # noqa: E402
from notebooks.models.utility import load_and_validate_final_data  # noqa: E402
from notebooks.preprocessing import assign_cohort  # noqa: E402

RESULTS_DIR = HERE / "results"

# Scores are evaluated as they come, with no sign flipping: each tool's intended direction is
# its own convention and guessing it would silently invent a result. A negative exp_med means
# the column ranks the opposite way to inhibition, and its magnitude is still the signal it
# carries. Only oligo_ai_score is a direct efficacy prediction; the rest are feature-like.
SANITY_FLOOR = 0.05

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--columns", nargs="+", default=None,
                    help="subset of the COMPETITION columns to score (default: all present)")
    args = ap.parse_args()

    # load_competition=True fetches and md5-checks the competition shards through the repo's
    # own cache, so the scores are read rather than recomputed.
    df, _ = load_and_validate_final_data(version="oligo", load_competition=True)
    _, test = common.split(assign_cohort(df))

    wanted = args.columns or COMPETITION
    payload, table_rows = [], {}
    for column in wanted:
        if column not in test.columns:
            print(f"{column}: not in the fetched shards, skipping")
            continue
        present = test[test[column].notna()]
        if present.empty:
            print(f"{column}: no scored test rows, skipping")
            continue
        if present[column].nunique(dropna=True) <= 1:
            # Same degeneracy check regen_common.validate_features makes: a constant column
            # cannot rank anything, and every per-experiment Spearman comes back undefined.
            print(f"{column}: degenerate (all NaN/constant) over the test rows, skipping")
            continue
        result = common.metrics_from(present[column].to_numpy(np.float64), present)
        suspect = abs(result["exp_med"]) < SANITY_FLOOR
        print(f"{column:22s} exp_med {result['exp_med']:+.3f}  on {len(present)} rows"
              + ("   <- no usable ranking signal" if suspect else ""))
        payload.append({"model": column, "kind": "published competitor",
                        "test_rows": len(present), "suspect": bool(suspect),
                        **{k: float(v) for k, v in result.items()}})
        table_rows[column] = result

    if not payload:
        raise SystemExit("no competitor columns scored")

    text = "\n".join([f"Published competitors on the held-out test of the frozen OligoAI split",
                      "", *common.metric_table("HELD-OUT TEST", table_rows)])
    print("\n" + text)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    (RESULTS_DIR / "competitors_reference.txt").write_text(text + "\n")
    (RESULTS_DIR / "competitors_reference.json").write_text(json.dumps(payload, indent=2))
    print(f"\nsaved -> {RESULTS_DIR / 'competitors_reference.txt'} / .json")


if __name__ == "__main__":
    main()
