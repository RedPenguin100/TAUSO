"""Score the deployed booster on the held-out test split and save the report.

The booster that `tauso setup-model` installs is the one the package scores with, so its
numbers belong in the repository next to the parameters it was trained under. Reports the
metric suite in notebooks.models.evaluate, all of it ranked against actual inhibition.

    python notebooks/models/evaluate_deployed.py [--version v1]
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from notebooks.models import common

from tauso.inference.scoring import DEFAULT_VERSION, MODEL_FILES, load_model


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--version", default=DEFAULT_VERSION, help="deployed model version to score")
    args = ap.parse_args()

    booster, features = load_model(args.version)
    df, _ = common.load_dataset()
    _, test = common.split(df)
    scores = common.metrics_from(common.predict(booster, test, features), test)

    lines = [
        f"Deployed model {args.version}: {MODEL_FILES[args.version]['filename']}",
        f"{len(features)} features, held-out test split, n={len(test)}",
        "",
        *common.metric_table("TEST", {args.version: scores}),
    ]
    common.write_report(lines, f"tauso_score_{args.version}_test_metrics.txt")


if __name__ == "__main__":
    main()
