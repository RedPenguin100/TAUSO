"""Step 1: download the frozen raw OligoAI flank-50 dataset from Zenodo (record 20794660)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from notebooks.consts import ORIGINAL_OLIGO_CSV_RAW
from tauso.cli_utils import download_with_progress, file_matches_hash, verify_hash_or_exit

ZENODO_RECORD = "20794660"
ZENODO_URL = (
    f"https://zenodo.org/api/records/{ZENODO_RECORD}/files/"
    "aso_inhibitions_21_08_25_incl_context_w_flank_50_df.csv.gz/content"
)
# md5 published by the Zenodo record; pins the exact frozen dataset.
EXPECTED_MD5 = "0efae5ae62d8d2c222c9caf5e01ced07"


def main():
    ORIGINAL_OLIGO_CSV_RAW.parent.mkdir(parents=True, exist_ok=True)

    if file_matches_hash(str(ORIGINAL_OLIGO_CSV_RAW), EXPECTED_MD5, algo="md5"):
        print(f"Raw flank-50 already present and verified -> {ORIGINAL_OLIGO_CSV_RAW}")
        return

    print(f"Downloading raw flank-50 from Zenodo -> {ORIGINAL_OLIGO_CSV_RAW}")
    download_with_progress(ZENODO_URL, str(ORIGINAL_OLIGO_CSV_RAW), label="Downloading raw flank-50")
    verify_hash_or_exit(str(ORIGINAL_OLIGO_CSV_RAW), EXPECTED_MD5, algo="md5")


if __name__ == "__main__":
    main()
