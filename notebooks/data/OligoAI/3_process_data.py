"""Step 3: compute binding-site structure, filter, and average technical replicates.

Produces OLIGO_CSV_PROCESSED_AVERAGED -- the label/split/metadata table the model reads.
Filtering (chemistry, missing gene, cell line, unmapped, ...) lives in
notebooks.preprocessing.process_oligo_data.
"""
import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from notebooks.consts import OLIGO_CSV_INDEXED, OLIGO_CSV_PROCESSED, OLIGO_CSV_PROCESSED_AVERAGED
from notebooks.preprocessing import process_oligo_data, process_oligo_data_rename
from tauso.data.consts import INHIBITION_PERCENT
from tauso.populate.calculators.calculator import Calculator


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cpus", type=int, default=1, help="CPUs for structure calculation")
    args = p.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s | %(message)s")

    data = process_oligo_data_rename(pd.read_csv(OLIGO_CSV_INDEXED))

    calc = Calculator(data=data, data_version="oligo", overwrite=True, cpus=args.cpus)
    calc.calculate_structure()

    processed = process_oligo_data(calc.data, strict_gapmer_patterns=False)
    processed.to_csv(OLIGO_CSV_PROCESSED)

    feature_cols = [c for c in processed.columns if c not in [INHIBITION_PERCENT, "index_oligo"]]
    averaged = processed.groupby(feature_cols, as_index=False, dropna=False).agg(
        {INHIBITION_PERCENT: "median", "index_oligo": "first"}
    )[["index_oligo"] + feature_cols + [INHIBITION_PERCENT]]
    averaged.to_csv(OLIGO_CSV_PROCESSED_AVERAGED, index=False)

    print(f"processed {len(processed):,} -> averaged {len(averaged):,} -> {OLIGO_CSV_PROCESSED_AVERAGED}")


if __name__ == "__main__":
    main()
