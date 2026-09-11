"""Fetch the DIAPH3 feature tilings from Zenodo, verifying each file's md5.

Three backbones over the same 498,327 candidates tiled across DIAPH3's pre-mRNA, scored for
OVCAR-8 / lipofection / 100 nM: a full phosphorothioate 5-10-5 2'-MOE gapmer and the two mixed
PS/PO backbones that cover 93.1% of the corpus. `_features.parquet` holds every feature the
model reads plus its score; `_ranked.csv` is the readable summary with the liability block.

  python notebooks/prediction/download_DIAPH3.py --ranked          # the summaries, 270 MB
  python notebooks/prediction/download_DIAPH3.py --backbone ps6    # one backbone, both files
  python notebooks/prediction/download_DIAPH3.py                   # everything, 3.8 GB
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tauso.cli._download import _ensure_zenodo_content_file

RECORD = "22693356"  # Computed features for DIAPH3 for the TAUSO model
OUT = Path(__file__).resolve().parent / "output" / "DIAPH3_ovcar8_lipo"

# backbone -> {kind: (filename, md5)}. `full` is the all-PS backbone; ps5 and ps6 carry 5 and 6
# phosphodiester linkages, both in the wings.
FILES = {
    "full": {
        "features": ("DIAPH3_ovcar8_lipo_full_features.parquet", "5ef405d58ec3eac2f0752b93d4b227e3"),
        "ranked": ("DIAPH3_ovcar8_lipo_full_ranked.csv", "9fc87085009b650ee819b89902b9ea44"),
    },
    "ps5": {
        "features": ("DIAPH3_ovcar8_lipo_ps5_features.parquet", "99cfb716e700dfbe1066157eb22913ed"),
        "ranked": ("DIAPH3_ovcar8_lipo_ps5_ranked.csv", "56efd0d237056547f72c4835d81bbe44"),
    },
    "ps6": {
        "features": ("DIAPH3_ovcar8_lipo_ps6_features.parquet", "772d9d34c13b54bad7daedca7489a6c1"),
        "ranked": ("DIAPH3_ovcar8_lipo_ps6_ranked.csv", "53b128bf928a6c63e1756c6f67e37650"),
    },
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--backbone", choices=sorted(FILES), action="append", help="repeatable; default all three")
    ap.add_argument("--ranked", action="store_true", help="only the summaries, not the feature tilings")
    ap.add_argument("--out", type=Path, default=OUT, help=f"destination directory (default {OUT})")
    ap.add_argument("--force", action="store_true", help="re-download even if the file is already there")
    args = ap.parse_args()

    kinds = ["ranked"] if args.ranked else ["features", "ranked"]
    args.out.mkdir(parents=True, exist_ok=True)
    for backbone in args.backbone or sorted(FILES):
        for kind in kinds:
            filename, md5 = FILES[backbone][kind]
            _ensure_zenodo_content_file(RECORD, filename, str(args.out / filename), md5, "md5", args.force)
    print(f"\n  -> {args.out}")


if __name__ == "__main__":
    main()
