"""Write tauso_score_<version>.finite.txt: the model features that are never missing in training.

Part of the feature set is conditionally defined -- codon scores exist only for windows inside a CDS,
cEt hybridization terms only for cEt chemistries, splice-junction distances only near a junction --
and the booster learned a missing-value branch for each. For the features listed here it never saw a
missing value, so a NaN in one of them at inference time is a failed computation rather than a real
input.

    python notebooks/models/write_finite_features.py [--version v1] [--run oligo]
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tauso.inference.predict import MODEL_DIR
from tauso.populate.feature_cache import cache_path_if_present


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--version", default="v1", help="model version to write the list for")
    ap.add_argument("--run", default="oligo", help="feature-store run to measure against")
    args = ap.parse_args()

    store = cache_path_if_present(args.run)
    if store is None:
        sys.exit(f"No feature store for run {args.run!r}; run `tauso setup-features` first.")

    features = (MODEL_DIR / f"tauso_score_{args.version}.features.txt").read_text().split()
    df = pd.read_parquet(store, columns=features)
    finite = sorted(f for f in features if not df[f].isna().any())

    out = MODEL_DIR / f"tauso_score_{args.version}.finite.txt"
    out.write_text("\n".join(finite) + "\n")
    print(f"{len(finite)} of {len(features)} features are never missing in {store} -> {out}")


if __name__ == "__main__":
    main()
