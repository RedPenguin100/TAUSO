"""Regenerate the OligoWalk competition features (OW_*) for the full-raw index, keyed by index_oligo."""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from notebooks.consts import OLIGO_CSV_INDEXED
from notebooks.notebook_utils import get_unique_genes
from notebooks.preprocessing import process_oligo_data_rename
from notebooks.features.feature_extraction import save_feature
from tauso.data.consts import CANONICAL_GENE_NAME
from tauso.genome.read_human_genome import get_locus_to_data_dict
from tauso.features.competition.oligowalk_utils import populate_oligowalk


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpus", type=int, default=None)
    args = ap.parse_args()

    data = process_oligo_data_rename(pd.read_csv(OLIGO_CSV_INDEXED))
    data = data[data[CANONICAL_GENE_NAME].notna()].reset_index(drop=True)
    gene_to_data = get_locus_to_data_dict(include_introns=True, gene_subset=get_unique_genes(data))

    data, features = populate_oligowalk(data, gene_to_data, num_workers=args.cpus)
    for f in features:
        if f != "error":
            save_feature(data, feature_name=f, overwrite=True, version="oligo")
    print("oligowalk done:", [f for f in features if f != "error"])


if __name__ == "__main__":
    main()
