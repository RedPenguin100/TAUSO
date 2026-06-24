"""Regenerate the sFold competition feature (sfold_accessibility) for the full-raw index."""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from notebooks.competitors.regen_common import load_gene_mapped, validate_and_report
from notebooks.competitors.sfold_utils import calculate_sfold_accessibility
from notebooks.features.feature_extraction import save_feature


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpus", type=int, default=None)
    args = ap.parse_args()

    data, gene_to_data = load_gene_mapped()
    data = calculate_sfold_accessibility(df=data, gene_to_data=gene_to_data, n_cores=args.cpus)
    validate_and_report(data, ["sfold_accessibility"])
    save_feature(data, "sfold_accessibility", overwrite=True, version="oligo")
    print("sfold done")


if __name__ == "__main__":
    main()
