"""
Calculate features for a dataset and save them to disk.

Each feature is written as an individual CSV. If the script crashes and is
re-run, already-written features are skipped automatically.

Usage (run from project root):
    python -m notebooks.features.calculate_features --dataset oligo --cpus 48
    python -m notebooks.features.calculate_features --dataset oligo --step hybridization
    python -m notebooks.features.calculate_features --dataset oligo --input /path/to/data.csv.gz
    python -m notebooks.features.calculate_features --dataset oligo --overwrite

Partitioned runs (split work across N SLURM jobs, then merge):
    python -m notebooks.features.calculate_features --dataset oligo --cpus 64 --partition 0 4
    python -m notebooks.features.calculate_features --dataset oligo --cpus 64 --partition 1 4
    python -m notebooks.features.calculate_features --dataset oligo --cpus 64 --partition 2 4
    python -m notebooks.features.calculate_features --dataset oligo --cpus 64 --partition 3 4
    python -m notebooks.features.calculate_features --dataset oligo --merge 4

SLURM note: run with `python -u` or set PYTHONUNBUFFERED=1 for live log output.
"""
import argparse
import logging
import shutil
import sys
from pathlib import Path

import pandas as pd

from notebooks.consts import OLIGO_CSV_INDEXED
from notebooks.data.OligoAI.utility import standardize_oligo_ai_data
from notebooks.features.feature_extraction import _get_saved_features_dir
from tauso.data.consts import CANONICAL_GENE
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


def partition_data(data: pd.DataFrame, k: int, n: int) -> pd.DataFrame:
    """Return the subset of data whose genes fall in partition k of n."""
    genes = sorted(data[CANONICAL_GENE].dropna().unique())
    partition_genes = set(genes[k::n])
    logger.info("Partition %d/%d: %d genes (of %d total)", k, n, len(partition_genes), len(genes))
    return data[data[CANONICAL_GENE].isin(partition_genes)].copy()


def get_partition_dir(main_feature_dir: Path, k: int, n: int) -> Path:
    """Return the feature directory for partition k of n."""
    return main_feature_dir.parent / f"{main_feature_dir.name}_p{k}of{n}"


def merge_partitions(main_feature_dir: Path, n_partitions: int, index_col: str) -> None:
    """Merge all partition feature directories into main_feature_dir."""
    main_feature_dir.mkdir(parents=True, exist_ok=True)
    for k in range(n_partitions):
        part_dir = get_partition_dir(main_feature_dir, k, n_partitions)
        if not part_dir.exists():
            logger.warning("Partition dir %s not found, skipping.", part_dir)
            continue
        for csv_file in sorted(part_dir.glob("*.csv")):
            dest = main_feature_dir / csv_file.name
            if not dest.exists():
                shutil.copy2(csv_file, dest)
                logger.info("Merged %s from partition %d.", csv_file.stem, k)
            else:
                existing_idx = pd.read_csv(dest, usecols=[index_col])[index_col]
                new_df = pd.read_csv(csv_file)
                new_rows = new_df[~new_df[index_col].isin(existing_idx)]
                if not new_rows.empty:
                    new_rows.to_csv(dest, mode='a', header=False, index=False)
                    logger.info("Appended %d rows for %s from partition %d.", len(new_rows), csv_file.stem, k)
    logger.info("Merge complete into %s.", main_feature_dir)


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
    parser.add_argument('--partition', nargs=2, type=int, metavar=('K', 'N'),
                        help='Run gene partition K of N (0-indexed). Saves to a partition-specific dir.')
    parser.add_argument('--merge', type=int, metavar='N',
                        help='Merge N partition directories into the main feature directory.')
    parser.add_argument('--log-level', default='INFO',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        help='Logging verbosity (default: INFO). Use DEBUG for full output.')
    args = parser.parse_args()

    if args.partition and args.merge:
        parser.error('--partition and --merge are mutually exclusive.')

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format='%(asctime)s %(levelname)s %(message)s',
        stream=sys.stdout,
        force=True,
    )
    logging.getLogger("numba").setLevel(logging.WARNING)

    config = DATASETS[args.dataset]
    main_feature_dir = _get_saved_features_dir(config['version'])
    index_col = f"index_{config['version']}"

    if args.merge:
        logger.info("Merging %d partitions into %s ...", args.merge, main_feature_dir)
        merge_partitions(main_feature_dir, n_partitions=args.merge, index_col=index_col)
        return

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

    if args.partition:
        k, n = args.partition
        data = partition_data(data, k=k, n=n)
        part_dir = get_partition_dir(main_feature_dir, k=k, n=n)
        get_feature_dir = lambda v, d=part_dir: d
        logger.info("Saving features to partition dir: %s", part_dir)
    else:
        get_feature_dir = _get_saved_features_dir

    calculator = Calculator(
        data=data,
        data_version=config['version'],
        overwrite=args.overwrite,
        cpus=args.cpus,
        get_feature_dir=get_feature_dir,
    )

    for step in steps:
        logger.info("Running %s ...", step)
        getattr(calculator, step)()
        logger.debug("Finished %s.", step)

    logger.info("Done.")


if __name__ == '__main__':
    main()
