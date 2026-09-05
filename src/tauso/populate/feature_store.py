"""Where computed features live on disk, and how they are found again.

A feature is stored either as a loose per-feature shard under ``<dir>/_patches/`` (or, for
the migration window, in ``<dir>`` itself) or as a column of the wide parquet cache. Shards
win, so a recomputed feature overrides the packed one.
"""

import logging
import os

import pandas as pd
import pyarrow.parquet as pq

from .feature_cache import FEATURE_CACHE_FILES, cache_path_if_present, loose_shard_dir, save_feature_internal

logger = logging.getLogger(__name__)

SHARD = "shard"
CACHE = "cache"


def index_column(version):
    """The index column for `version`, or None when features are not being stored."""
    if version is None:
        return None
    if version in FEATURE_CACHE_FILES:
        return FEATURE_CACHE_FILES[version]["index_col"]
    return f"index_{version}"


class FeatureStore:
    """Disk-backed feature storage for one data version.

    `overwrite` makes every feature report as missing, so a run recomputes and replaces
    what is already stored.
    """

    def __init__(self, version, get_feature_dir=None, overwrite=False):
        self.version = version
        self.overwrite = overwrite
        self.index = index_column(version)
        self._get_feature_dir = get_feature_dir
        self._wide = None

    @property
    def enabled(self):
        """Whether features can be read from or written to disk at all."""
        return self._get_feature_dir is not None and self.version is not None

    @property
    def base_dir(self):
        if not self.enabled:
            raise ValueError("no feature directory: the store was built without get_feature_dir")
        return self._get_feature_dir(self.version)

    def _wide_cache(self):
        """`(column names, path)` of the wide cache, or `(set(), None)`.

        Read once: the wide cache is never written during a run, unlike the shards.
        """
        if self._wide is None:
            path = cache_path_if_present(self.version) if self.version else None
            self._wide = (set(pq.read_schema(path).names), path) if path else (set(), None)
        return self._wide

    def locate(self, feature):
        """`(SHARD, path)`, `(CACHE, path)` or None. Shards win.

        Shard paths are stat-ed on every call because a step writes them as it runs.
        """
        if not self.enabled:
            return None
        base = self.base_dir
        for parent in (loose_shard_dir(base), base):
            for ext in (".parquet", ".csv"):
                path = os.path.join(parent, f"{feature}{ext}")
                if os.path.exists(path):
                    return (SHARD, path)
        cache_cols, cache = self._wide_cache()
        if feature in cache_cols:
            return (CACHE, cache)
        return None

    def missing(self, features):
        """Those of `features` that are not stored, in the order given."""
        if self.overwrite or not self.enabled:
            return list(features)
        missing = []
        for feature in features:
            if self.locate(feature) is None:
                missing.append(feature)
            else:
                logger.debug("Skipping '%s': already present.", feature)
        return missing

    def read(self, features, index_values=None):
        """One frame holding `index` and `features`, read with a single pass per file.

        Every feature in the wide cache is fetched in one read, rather than reopening the
        parquet once per column. `index_values` anchors the result to the rows the caller
        wants: a source covering other rows would otherwise contribute them, and the columns
        it does not cover would widen to float to hold their NaN.
        """
        if self.index is None:
            raise ValueError(f"cannot load {list(features)} from disk without a data version")
        if not self.enabled:
            raise ValueError(f"cannot load {list(features)} without a feature directory")

        by_source = {}
        for feature in dict.fromkeys(features):  # a repeat would be read, and merged, twice
            source = self.locate(feature)
            if source is None:
                raise FileNotFoundError(
                    f"No shard or cache column for '{feature}' in {self.base_dir} (.parquet, .csv, or wide cache)."
                )
            by_source.setdefault(source, []).append(feature)

        frame = None if index_values is None else pd.DataFrame({self.index: index_values})
        for (kind, path), columns in by_source.items():
            if kind == CACHE:
                part = pd.read_parquet(path, columns=[self.index, *columns])
            elif path.endswith(".parquet"):
                part = pd.read_parquet(path)[[self.index, *columns]]
            else:
                part = pd.read_csv(path)[[self.index, *columns]]
            frame = part if frame is None else frame.merge(part, on=self.index, how="left")
            logger.debug("Read %s from %s.", columns, kind)
        return frame

    def save(self, data, feature):
        """Write `feature` as a shard. A store without a feature directory writes nothing."""
        if self._get_feature_dir is None:
            return
        save_feature_internal(
            data,
            feature_name=feature,
            overwrite=self.overwrite,
            version=self.version,
            saved_dir_func=self._get_feature_dir,
        )
