"""Step 2: assign the canonical gene from the patent's target_gene, with curated corrections.

The gene comes from the patent label, never re-derived from alignment: an ASO whose patent targets
gene X must not be silently reassigned to a gene Y it happens to align to. Rows with no target_gene,
or whose gene is an unrecoverable correction (mapped to NA), are left NaN and dropped in
3_process_data. (See curate_gene_labels.py for how the correction map was derived.)
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from notebooks.consts import ORIGINAL_OLIGO_CSV, ORIGINAL_OLIGO_CSV_WITH_CANONICAL
from notebooks.data.OligoAI.gene_corrections import MANUAL_CANONICAL_MAPPING
from tauso.data.consts import CANONICAL_GENE_NAME


def main():
    df = pd.read_csv(ORIGINAL_OLIGO_CSV)
    df[CANONICAL_GENE_NAME] = df["target_gene"].replace(MANUAL_CANONICAL_MAPPING)
    n_nan = int(df[CANONICAL_GENE_NAME].isna().sum())
    df.to_csv(ORIGINAL_OLIGO_CSV_WITH_CANONICAL, index=False)
    print(f"assigned canonical gene for {len(df):,} rows "
          f"({n_nan:,} left NaN -> dropped in 3_process_data)")


if __name__ == "__main__":
    main()
