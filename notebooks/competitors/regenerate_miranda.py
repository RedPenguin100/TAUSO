"""Regenerate the miRanda competition features (miranda_score/energy) for the full-raw index."""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from notebooks.competitors.regen_common import load_gene_mapped, validate_and_report
from notebooks.competitors.miranda_utils import get_populated_df_with_miranda_features
from notebooks.features.feature_extraction import save_feature


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpus", type=int, default=None)
    args = ap.parse_args()

    data, gene_to_data = load_gene_mapped()
    data = get_populated_df_with_miranda_features(data, gene_to_data, max_workers=args.cpus)
    validate_and_report(data, ["miranda_score", "miranda_energy"])
    save_feature(data, "miranda_score", overwrite=True, version="oligo")
    save_feature(data, "miranda_energy", overwrite=True, version="oligo")
    print("miranda done")


if __name__ == "__main__":
    main()
