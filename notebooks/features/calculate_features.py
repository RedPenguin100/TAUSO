"""
Calculate features for a dataset and save them to disk.

Each feature is written as an individual CSV. If the script crashes and is
re-run, already-written features are skipped automatically.

Usage (run from project root):
    python -m notebooks.features.calculate_features --dataset oligo --cpus 48
    python -m notebooks.features.calculate_features --dataset oligo --step hybridization
    python -m notebooks.features.calculate_features --dataset oligo --input /path/to/data.csv.gz
    python -m notebooks.features.calculate_features --dataset oligo --overwrite

SLURM note: run with `python -u` or set PYTHONUNBUFFERED=1 for live log output.
"""
import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

from notebooks.consts import OLIGO_CSV_INDEXED
from notebooks.data.OligoAI.utility import standardize_oligo_ai_data
from notebooks.features.feature_extraction import _get_saved_features_dir
from tauso.populate.calculators.calculator import Calculator

logger = logging.getLogger(__name__)


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
    parser.add_argument('--step', default=None,
                        help='Run only one step (e.g. hybridization or calculate_hybridization). '
                             'Runs all steps if omitted.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        help='Logging verbosity (default: INFO). Use DEBUG for full output.')
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format='%(asctime)s %(levelname)s %(message)s',
        stream=sys.stdout,
        force=True,
    )

    config = DATASETS[args.dataset]
    steps = config['steps']

    if args.step:
        canonical = args.step if args.step.startswith('calculate_') else f'calculate_{args.step}'
        if canonical not in steps:
            valid = ', '.join(s.removeprefix('calculate_') for s in steps)
            parser.error(f"Unknown step '{args.step}'. Valid steps for '{args.dataset}': {valid}")
        steps = [canonical]

    csv_path = args.input or config['default_input']

    logger.info("Loading %s data from %s ...", args.dataset, csv_path)
    data = config['loader'](csv_path)
    logger.info("Loaded %d rows.", len(data))

    calculator = Calculator(
        data=data,
        data_version=config['version'],
        overwrite=args.overwrite,
        cpus=args.cpus,
        get_feature_dir=_get_saved_features_dir,
    )

    for step in steps:
        logger.info("Running %s ...", step)
        getattr(calculator, step)()
        logger.debug("Finished %s.", step)

    logger.info("Done.")


if __name__ == '__main__':
    main()
