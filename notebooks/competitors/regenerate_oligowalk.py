"""Regenerate the OligoWalk competition features (OW_*) for the full-raw index, keyed by index_oligo."""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from notebooks.competitors.regen_common import load_gene_mapped, validate_and_report
from notebooks.features.feature_extraction import save_feature
from tauso.features.competition.oligowalk_utils import populate_oligowalk


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpus", type=int, default=None)
    args = ap.parse_args()

    data, gene_to_data = load_gene_mapped()
    data, features = populate_oligowalk(data, gene_to_data, num_workers=args.cpus)
    feats = [f for f in features if f != "error"]
    validate_and_report(data, feats)
    for f in feats:
        save_feature(data, feature_name=f, overwrite=True, version="oligo")
    print("oligowalk done:", feats)


if __name__ == "__main__":
    main()
