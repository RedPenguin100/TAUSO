"""Step 2.5: add a stable sequential index_oligo key, so every downstream merge is order-independent."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from notebooks.consts import OLIGO_CSV_INDEXED, ORIGINAL_OLIGO_CSV_WITH_CANONICAL


def main():
    df = pd.read_csv(ORIGINAL_OLIGO_CSV_WITH_CANONICAL)
    df["index_oligo"] = range(1, len(df) + 1)
    df.to_csv(OLIGO_CSV_INDEXED, index=False)
    print(f"indexed {len(df):,} rows -> {OLIGO_CSV_INDEXED}")


if __name__ == "__main__":
    main()
