"""Step 1: download the frozen raw OligoAI flank-50 dataset from Zenodo (record 20794660)."""
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from notebooks.consts import ORIGINAL_OLIGO_CSV_RAW

ZENODO_URL = (
    "https://zenodo.org/api/records/20794660/files/"
    "aso_inhibitions_21_08_25_incl_context_w_flank_50_df.csv.gz/content"
)


def main():
    ORIGINAL_OLIGO_CSV_RAW.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading raw flank-50 from Zenodo -> {ORIGINAL_OLIGO_CSV_RAW}")
    urllib.request.urlretrieve(ZENODO_URL, ORIGINAL_OLIGO_CSV_RAW)


if __name__ == "__main__":
    main()
