"""Regenerate the ClinASO competition feature (ClinASO_score) for the full-raw index."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from notebooks.competitors.ClinASO.clinaso_utils import populate_clinaso
from notebooks.competitors.regen_common import load_indexed, validate_features
from notebooks.features.feature_extraction import save_feature


def main():
    data = load_indexed()
    data, features = populate_clinaso(data)
    validate_features(data, features)
    for feature in features:
        save_feature(data, feature, overwrite=True, version="oligo")
    print("clinaso done")


if __name__ == "__main__":
    main()
