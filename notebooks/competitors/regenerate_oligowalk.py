"""Regenerate the OligoWalk competition features (OW_*) for the full-raw index, keyed by index_oligo."""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from notebooks.competitors.regen_common import load_gene_mapped, validate_features
from notebooks.features.feature_extraction import save_feature
from tauso.features.competition.oligowalk_utils import (
    FLAGS_DEFAULT,
    FLAGS_DNA_TARGET_STRUCTURE,
    populate_oligowalk,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpus", type=int, default=None)
    ap.add_argument(
        "--target-structure",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="DNA oligomer parameters and a target refolded after binding (OligoWalk -d -m 2 -s 2, about 40 s per ASO); "
        "--no-target-structure runs -m 1 with RNA oligomer parameters and the target unfolded.",
    )
    args = ap.parse_args()

    data, gene_to_data = load_gene_mapped()
    flags = FLAGS_DNA_TARGET_STRUCTURE if args.target_structure else FLAGS_DEFAULT
    data, features = populate_oligowalk(data, gene_to_data, num_workers=args.cpus, flags=flags)
    feats = [f for f in features if f != "error"]
    # OW_Break_Target is constant 0 under OligoWalk -m 1; skip it in the degeneracy guard (still saved).
    guarded = feats if args.target_structure else [f for f in feats if f != "OW_Break_Target"]
    validate_features(data, guarded)
    for f in feats:
        save_feature(data, feature_name=f, overwrite=True, version="oligo")
    print("oligowalk done:", feats)


if __name__ == "__main__":
    main()
