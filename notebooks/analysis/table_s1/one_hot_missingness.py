"""Table S1: missingness of the terminal one-hot columns, and the rule that produces it.

A terminal one-hot position is left NaN when it falls outside its own half of the ASO, so
short designs carry no value at the inner positions. This reports the count per column and
checks it against the encoder's midpoint test.
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from notebooks.models.common import load_dataset
from tauso.data.consts import ASO_SEQUENCE
from tauso.populate.populate_sequence import TERMINAL_OHE_N


def min_length(prefix: str, position: int) -> int:
    """Shortest ASO for which ``prefix``/``position`` lies in its own half of the sequence.

    Mirrors the encoder: the 5' window keeps position i while ``2 * i <= n - 1`` and the 3'
    window while ``2 * (n - 1 - i) > n - 1``, so the 3' threshold is one nucleotide higher.
    """
    return 2 * position + (1 if prefix == "ohe_pos" else 2)


def missingness(df: pd.DataFrame) -> pd.DataFrame:
    """One row per terminal one-hot column: observed NaN count and the count the rule predicts."""
    lengths = df[ASO_SEQUENCE].str.len()
    rows = []
    for prefix in ("ohe_pos", "ohe_3p"):
        for position in range(TERMINAL_OHE_N):
            needs = min_length(prefix, position)
            for base in "ACGT":
                column = f"{prefix}{position}_{base}"
                rows.append(
                    {
                        "column": column,
                        "position": f"{prefix}{position}_*",
                        "min_aso_length": needs,
                        "observed": int(df[column].isna().sum()),
                        "predicted": int((lengths < needs).sum()),
                    }
                )
    out = pd.DataFrame(rows)
    out["pct"] = (100 * out["observed"] / len(df)).round(2)
    return out


def main():
    df, _ = load_dataset()
    per_column = missingness(df)

    disagree = per_column[per_column["observed"] != per_column["predicted"]]
    if not disagree.empty:
        raise AssertionError(f"midpoint rule does not explain the NaN mask:\n{disagree}")

    print(f"{len(df):,} rows, {len(per_column)} terminal one-hot columns")
    print(f"ASO length: {df[ASO_SEQUENCE].str.len().min()}-{df[ASO_SEQUENCE].str.len().max()} nt\n")
    grouped = per_column.groupby(
        ["position", "min_aso_length", "observed", "pct"], as_index=False
    ).agg(columns=("column", "count"))
    print(grouped.sort_values("min_aso_length").to_string(index=False))
    print("\nevery column agrees with the midpoint rule")


if __name__ == "__main__":
    main()
