"""
Calculate features for a dataset and save them to disk.

Each feature is written as an individual CSV. If the script crashes and is
re-run, already-written features are skipped automatically.

Usage:
    python notebooks/features/calculate_features.py --dataset oligo --cpus 48
    python notebooks/features/calculate_features.py --dataset oligo --input /path/to/data.csv.gz
    python notebooks/features/calculate_features.py --dataset oligo --overwrite
"""
import argparse
from pathlib import Path

import pandas as pd

from notebooks.consts import OLIGO_CSV_INDEXED
from notebooks.data.OligoAI.utility import standardize_oligo_ai_data
from tauso.populate.calculators.calculator import Calculator


def _load_oligo(csv_path: Path) -> pd.DataFrame:
    return standardize_oligo_ai_data(pd.read_csv(csv_path))


DATASETS = {
    'oligo': {
        'loader': _load_oligo,
        'version': 'oligo',
        'default_input': OLIGO_CSV_INDEXED,
        'steps': [
            'calculate_structure',
            'calculate_sense_accessibility',
            'calculate_sequence_one_hot',
            'calculate_sequence_chemistry',
            'calculate_modification',
            'calculate_hybridization',
            'calculate_backbone_features',
            'calculate_ribo_seq',
            'calculate_off_target_general',
            'calculate_off_target_single',
            'calculate_mrna_halflife',
            'calculate_rbp',
            'calculate_off_target_specific',
        ],
    },
    # 'asoptimizer': { ... }  # TODO: add when needed
}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--dataset', required=True, choices=list(DATASETS),
                        help='Which dataset to process')
    parser.add_argument('--input', type=Path, default=None,
                        help='Path to input CSV (overrides the dataset default)')
    parser.add_argument('--cpus', type=int, default=32,
                        help='CPUs for parallel steps (default: 32)')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite existing feature files')
    args = parser.parse_args()

    config = DATASETS[args.dataset]
    csv_path = args.input or config['default_input']

    print(f"Loading {args.dataset} data from {csv_path} ...")
    data = config['loader'](csv_path)
    print(f"Loaded {len(data)} rows.")

    calculator = Calculator(
        data=data,
        data_version=config['version'],
        overwrite=args.overwrite,
        cpus=args.cpus,
    )

    for step in config['steps']:
        print(f"\n[{step}]")
        getattr(calculator, step)()

    print("\nDone.")


if __name__ == '__main__':
    main()
