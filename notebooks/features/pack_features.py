"""Pack a CSV feature cache into the single wide Parquet distributed as the feature cache.

The ``notebooks/saved_features_<run>`` caches are large (the ``oligo`` run is ~1.7 GB of
single-column CSVs) and should not live in git. This tool outer-joins every feature into one
wide ``<run>_features.parquet`` (zstd, with native dtypes/nulls) for upload to Zenodo. At load
time it is just one (multi-column) parquet in the store folder alongside any per-feature dev
shards; ``tauso setup-features`` downloads it into ``TAUSO_DATA_DIR/features/<run>/``.

Release/dev tool, not part of the importable library -- it lives under notebooks/ and reuses
the per-feature dtype map from feature_extraction so the packed Parquet matches the CSVs.

Usage:
    python -m notebooks.features.pack_features --run oligo
    python -m notebooks.features.pack_features --source <dir> --run oligo --out build/
"""

import argparse
import concurrent.futures
import logging
from pathlib import Path

import pandas as pd

from notebooks.features.feature_extraction import _get_saved_features_dir, get_dtype_for_feature
from tauso.cli_utils import md5_file


def pack_features(source_dir: Path, out_dir: Path, run: str, workers: int, row_group_size: int = 20000) -> Path:
    """Outer-join every feature CSV in ``source_dir`` into one wide ``<run>_features.parquet``.

    ``row_group_size`` controls the parquet row-group granularity: smaller groups let
    ``iter_all_features`` stream the cache with peak memory ~= one row group.
    """
    index_col = f"index_{run}"
    csvs = sorted(p for p in source_dir.glob("*.csv") if not p.name.startswith("."))
    if not csvs:
        raise FileNotFoundError(f"No CSV features found in {source_dir}")

    def read_one(p: Path) -> pd.DataFrame:
        return pd.read_csv(
            p,
            index_col=index_col,
            dtype=get_dtype_for_feature(p.name, index_col),
            engine="pyarrow",
            dtype_backend="pyarrow",
        )

    logging.info("Reading %d feature CSVs ...", len(csvs))
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
        frames = list(ex.map(read_one, csvs))

    logging.info("Outer-joining into one wide table ...")
    wide = pd.concat(frames, axis=1, join="outer", copy=False).reset_index()

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{run}_features.parquet"
    logging.info("Writing %d columns x %d rows -> %s", wide.shape[1] - 1, len(wide), out_path.name)
    wide.to_parquet(out_path, compression="zstd", index=False, row_group_size=row_group_size)
    return out_path


def _dir_size(path: Path) -> int:
    return sum(f.stat().st_size for f in path.glob("*") if f.is_file())


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run", default="oligo", help="Feature-set run / dataset tag (default: oligo).")
    parser.add_argument(
        "--source", type=Path, default=None, help="Source CSV cache dir (default: the saved_features_<run> dir)."
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("feature_pack_build"),
        help="Build output dir for the parquet (default: ./feature_pack_build).",
    )
    parser.add_argument("--workers", type=int, default=8, help="Parallel CSV readers (default: 8).")
    parser.add_argument(
        "--row-group-size", type=int, default=20000, help="Parquet row-group size in rows (default: 20000)."
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S")

    source_dir = args.source or Path(_get_saved_features_dir(args.run))
    if not source_dir.is_dir():
        parser.error(f"Source dir does not exist: {source_dir}")

    out_path = pack_features(source_dir, args.out, args.run, args.workers, args.row_group_size)

    csv_bytes = _dir_size(source_dir)
    parquet_bytes = out_path.stat().st_size
    parquet_md5 = md5_file(str(out_path))

    logging.info("Done.")
    logging.info("  source CSVs : %8.1f MB  (%s)", csv_bytes / 1e6, source_dir)
    logging.info(
        "  parquet     : %8.1f MB  (%.1fx smaller)  %s",
        parquet_bytes / 1e6,
        csv_bytes / max(parquet_bytes, 1),
        out_path,
    )
    logging.info("  md5         : %s", parquet_md5)
    logging.info("")
    logging.info("After uploading to Zenodo, set in src/tauso/populate/feature_cache.py:")
    logging.info('    ZENODO_FEATURES_RECORD = "<record id>"')
    logging.info(
        '    FEATURE_CACHE_FILES["%s"] = {"filename": "%s", "md5": "%s"}', args.run, out_path.name, parquet_md5
    )


if __name__ == "__main__":
    main()
