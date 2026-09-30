"""The TAUSO row: the shipped booster scored on the same held-out test rows.

Loads the repo's dataset and frozen split, scores the test side with the deployed
`tauso_score_v1` booster, and evaluates it with the same `common.metrics_from` every other row
in the comparison uses. Nothing is trained here.

Run (needs the tauso package and TAUSO_DATA_DIR, i.e. the tauso_claude_ws environment):
  python notebooks/competitors/oligogym/score_tauso.py
"""

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
sys.path[:0] = [str(REPO_ROOT), str(REPO_ROOT / "src")]

from notebooks.models import common  # noqa: E402

from tauso.inference.scoring import MODEL_DIR, predict  # noqa: E402

RESULTS_DIR = HERE / "results"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--version", default="v1", help="model version in tauso.inference.scoring")
    args = ap.parse_args()

    features = (MODEL_DIR / f"tauso_score_{args.version}.features.txt").read_text().splitlines()
    df, _ = common.load_dataset()
    _, test = common.split(df)
    print(f"scoring {len(test)} test rows with {len(features)} model features", flush=True)

    result = common.metrics_from(predict(test, version=args.version, strict=False), test)

    label = f"tauso_score_{args.version}"
    lines = [f"TAUSO shipped model on the held-out test of the frozen OligoAI split"
             f"   test={len(test)}   {len(features)} features", "",
             *common.metric_table("HELD-OUT TEST (genome-derived features)", {label: result})]
    text = "\n".join(lines)
    print("\n" + text)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    (RESULTS_DIR / "tauso_reference.txt").write_text(text + "\n")
    (RESULTS_DIR / "tauso_reference.json").write_text(json.dumps(
        [{"model": label, "kind": "genome-derived", "n_features": len(features),
          "test_rows": len(test), **{k: float(v) for k, v in result.items()}}], indent=2))
    print(f"\nsaved -> {RESULTS_DIR / 'tauso_reference.txt'} / .json")


if __name__ == "__main__":
    main()
