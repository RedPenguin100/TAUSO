"""Step 2: assign the canonical gene from the patent's target_gene, with curated corrections.

The gene comes from the patent label, never re-derived from alignment: an ASO whose patent targets
gene X must not be silently reassigned to a gene Y it happens to align to. When target_gene is blank
but the patent still names the target in target_mrna, that label is used as a fallback, but only when
it is a curated alias in the correction map (so an unmapped free-text label can never become a gene).
A few aliases name different genes in different source tables of one patent; those are resolved by
custom_id via AMBIGUOUS_ALIAS_BY_CUSTOM_ID before the generic fallback.
Rows with no usable label, or whose gene is an unrecoverable correction (mapped to NA), are left NaN
and dropped in 3_process_data. (See curate_gene_labels.py for how the correction map was derived.)
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from notebooks.consts import ORIGINAL_OLIGO_CSV, ORIGINAL_OLIGO_CSV_WITH_CANONICAL
from notebooks.data.OligoAI.gene_corrections import AMBIGUOUS_ALIAS_BY_CUSTOM_ID, MANUAL_CANONICAL_MAPPING
from tauso.data.consts import CANONICAL_GENE_NAME


def main():
    df = pd.read_csv(ORIGINAL_OLIGO_CSV)
    n_ambig = 0
    if "custom_id" in df.columns:                                  # ambiguous aliases: gene depends on source table
        cid = df["custom_id"].astype(str)
        for token, gene in AMBIGUOUS_ALIAS_BY_CUSTOM_ID.items():
            m = df["target_gene"].isna() & cid.str.contains(token, regex=False, na=False)
            df.loc[m, "target_gene"] = gene
            n_ambig += int(m.sum())
    fallback = df["target_gene"].isna() & df["target_mrna"].isin(list(MANUAL_CANONICAL_MAPPING))
    df.loc[fallback, "target_gene"] = df.loc[fallback, "target_mrna"]
    df[CANONICAL_GENE_NAME] = df["target_gene"].replace(MANUAL_CANONICAL_MAPPING)
    n_nan = int(df[CANONICAL_GENE_NAME].isna().sum())
    df.to_csv(ORIGINAL_OLIGO_CSV_WITH_CANONICAL, index=False)
    print(f"assigned canonical gene for {len(df):,} rows "
          f"({int(fallback.sum()):,} recovered from target_mrna alias; "
          f"{n_ambig:,} from table-specific ambiguous alias; "
          f"{n_nan:,} left NaN -> dropped in 3_process_data)")


if __name__ == "__main__":
    main()
