"""Regenerate the PFRED competition features (PFRED_SVM/PFRED_PLS) for the full-raw index.

Pre-conditions (checked at runtime; each prints how to fix itself):
  1. PFREDIntegration submodule initialized:
       git submodule update --init notebooks/competitors/PFRED/PFREDIntegration
  2. the 'pfred' Docker container (image tauso/pfred:v1) running:
       docker start pfred                              # if it already exists (stopped)
       docker run -d -t --name pfred tauso/pfred:v1    # to create it the first time
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tauso.data.consts import ASO_SEQUENCE
from notebooks.competitors.PFRED.pfred_glue import populate_pfred, validate_docker_container
from notebooks.competitors.regen_common import load_indexed, validate_features
from notebooks.features.feature_extraction import save_feature


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--container", default="pfred")
    args = ap.parse_args()

    if not validate_docker_container(container_name=args.container, verbose=True):
        sys.exit("PFRED container not ready -- start it (see above) and re-run.")

    data = load_indexed()
    data, features = populate_pfred(data, seq_col=ASO_SEQUENCE, container_name=args.container)
    validate_features(data, list(features))
    for f in features:
        save_feature(data, f, overwrite=True, version="oligo")
    print("pfred done:", list(features))


if __name__ == "__main__":
    main()
