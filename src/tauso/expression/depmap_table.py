"""A DepMap expression table on disk, written and read without knowing how it is stored.

Every function takes the table's CSV path (the name DepMap publishes it under), whether or not
the CSV is still there. The CSVs are tens of thousands of columns wide, and a Parquet writer's
memory follows the width it is given, so a table is stored as a directory of narrow Parquet
parts: the first holds the profile columns, each later one COLUMNS_PER_PART measurement columns,
and every part the same rows in the same order. A data dir converted before the split holds one
Parquet file instead, which reads the same way.
"""

import hashlib
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

from ..cli_utils import sha256_file

DEPMAP_PROFILE_COLUMNS = frozenset(
    {"Unnamed: 0", "", "SequencingID", "ModelID", "IsDefaultEntryForModel", "ModelConditionID", "IsDefaultEntryForMC"}
)
"""The columns of a DepMap expression table that describe the sequencing run, not a measurement."""

COLUMNS_PER_PART = 2000
"""Measurement columns per Parquet part, and per read. A Parquet writer holds buffers for every
column until a row group is written, so narrow parts keep its memory small."""

CSV_BLOCK_BYTES = 1 << 28
"""Text read at once when pyarrow converts a part."""

ROWS_PER_GROUP = 600
"""Cell lines per Parquet row group when pyarrow converts a part."""


def _parts_dir(csv_path):
    path = Path(csv_path)
    return path.with_name(path.stem + "_parts")


def _single_file(csv_path):
    return Path(csv_path).with_suffix(".parquet")


def _location(csv_path):
    """Where the table is stored: its parts directory, or the single file of an older data dir."""
    single = _single_file(csv_path)
    return single if single.exists() and not _parts_dir(csv_path).is_dir() else _parts_dir(csv_path)


def _files(csv_path):
    """The table's Parquet files in column order; empty if it has not been written."""
    location = _location(csv_path)
    if location.is_dir():
        return sorted(str(p) for p in location.glob("part-*.parquet"))
    return [str(location)] if location.exists() else []


def _sidecar(csv_path):
    location = _location(csv_path)
    return location.with_name(location.name + ".sha256")


def _measurement_columns(schema):
    """All but the profile columns and a stored pandas index."""
    index_cols = {c for c in (schema.pandas_metadata or {}).get("index_columns", []) if isinstance(c, str)}
    return [c for c in schema.names if c not in DEPMAP_PROFILE_COLUMNS and c not in index_cols]


def table_exists(csv_path):
    """Whether the table has been written."""
    return bool(_files(csv_path))


def column_names(csv_path):
    """The table's measurement columns, in order."""
    return [c for path in _files(csv_path) for c in _measurement_columns(pq.read_schema(path))]


def _split_csv_by_columns(csv_path, out_dir):
    """Split the CSV into narrow CSVs of the same rows, and return their paths in column order.

    The split is plain text, one line at a time, so no parser ever sees the full width.
    """
    with open(csv_path) as source:
        header = source.readline().rstrip("\r\n").split(",")
        profile = [i for i, c in enumerate(header) if c in DEPMAP_PROFILE_COLUMNS]
        measurements = [i for i, c in enumerate(header) if c not in DEPMAP_PROFILE_COLUMNS]
        groups = [profile] + [
            measurements[start : start + COLUMNS_PER_PART] for start in range(0, len(measurements), COLUMNS_PER_PART)
        ]
        paths = [Path(out_dir) / f"part-{k:03d}.csv" for k in range(len(groups))]
        outs = [open(path, "w") for path in paths]
        try:
            for out, group in zip(outs, groups):
                out.write(",".join(header[i] for i in group) + "\n")
            for line in source:
                fields = line.rstrip("\r\n").split(",")
                for out, group in zip(outs, groups):
                    out.write(",".join([fields[i] for i in group]) + "\n")
        finally:
            for out in outs:
                out.close()
    return paths


def _pyarrow_csv_to_parquet(csv_path, parquet_path):
    with open(csv_path) as handle:
        header = handle.readline().rstrip("\n").split(",")
    # Typing the measurement columns up front stops a block of blanks being read as text,
    # which would clash with the float blocks either side of it.
    reader = pacsv.open_csv(
        csv_path,
        read_options=pacsv.ReadOptions(block_size=CSV_BLOCK_BYTES),
        convert_options=pacsv.ConvertOptions(
            column_types={c: pa.float64() for c in header if c not in DEPMAP_PROFILE_COLUMNS}
        ),
    )
    pending = []
    with pq.ParquetWriter(parquet_path, reader.schema) as writer:
        for batch in reader:
            pending.append(batch)
            if sum(b.num_rows for b in pending) >= ROWS_PER_GROUP:
                writer.write_table(pa.Table.from_batches(pending))
                pending = []
        if pending:
            writer.write_table(pa.Table.from_batches(pending))


def _pandas_csv_to_parquet(csv_path, parquet_path):
    pd.read_csv(csv_path).to_parquet(parquet_path, index=False)


PARSERS = {"pyarrow": _pyarrow_csv_to_parquet, "pandas": _pandas_csv_to_parquet}


def remove_table(csv_path):
    """Remove the table and its recorded hash, however it is stored."""
    for location in (_parts_dir(csv_path), _single_file(csv_path)):
        if location.is_dir():
            shutil.rmtree(location)
        else:
            location.unlink(missing_ok=True)
        location.with_name(location.name + ".sha256").unlink(missing_ok=True)


def write_table(csv_path, parser):
    """Convert the CSV to the table, replacing any earlier one, record its hash, and remove the CSV.

    parser is "pyarrow" or "pandas". The two round the last bit of some values differently, so
    each table is parsed by the one its outputs have always been built from.
    """
    convert = PARSERS[parser]
    directory = _parts_dir(csv_path)
    partial = directory.with_name(directory.name + ".partial")
    shutil.rmtree(partial, ignore_errors=True)
    partial.mkdir()
    for part_csv in _split_csv_by_columns(csv_path, partial):
        convert(str(part_csv), str(part_csv.with_suffix(".parquet")))
        part_csv.unlink()
    remove_table(csv_path)
    # Named only once it is whole, so a killed conversion is not mistaken for a finished one.
    os.replace(partial, directory)
    os.remove(csv_path)
    record_sha256(csv_path)


def table_sha256(csv_path):
    """The table's SHA-256: of its single file, or of its parts' names and hashes."""
    location = _location(csv_path)
    if not location.is_dir():
        return sha256_file(str(location))
    listing = "".join(f"{Path(p).name} {sha256_file(p)}\n" for p in _files(csv_path))
    return hashlib.sha256(listing.encode()).hexdigest()


def recorded_sha256(csv_path):
    """The hash recorded when the table was written, or None if none was."""
    sidecar = _sidecar(csv_path)
    return sidecar.read_text().strip() if sidecar.exists() else None


def record_sha256(csv_path):
    """Record the table's hash, so a later run can tell it has changed on disk."""
    _sidecar(csv_path).write_text(table_sha256(csv_path))


def read_rows(csv_path, model_ids):
    """The wanted cell lines' measurements, as (model ids found, measurement columns, values).

    A model can carry several sequencing profiles; the one DepMap flags as the default is read,
    and a model with two left is read from its first. Missing values read as 0. Only the wanted
    rows are kept, a batch of columns at a time, so what stays resident is the cohort, not the
    table.
    """
    paths = _files(csv_path)
    first_part = pq.ParquetFile(paths[0])
    names = first_part.schema.names
    meta = first_part.read(columns=[c for c in names if c in DEPMAP_PROFILE_COLUMNS]).to_pandas()
    model_col = "ModelID" if "ModelID" in meta.columns else meta.columns[0]
    wanted = meta[model_col].isin(model_ids)
    # DepMap flags the default profile with the strings "Yes"/"No".
    if "IsDefaultEntryForModel" in meta.columns:
        wanted &= meta["IsDefaultEntryForModel"] == "Yes"
    rows = np.flatnonzero(wanted.to_numpy())
    _, first = np.unique(meta[model_col].to_numpy()[rows], return_index=True)
    rows = rows[np.sort(first)]

    # One part is open at a time: every open part holds its footer in memory.
    part_cols = [_measurement_columns(first_part.schema_arrow)]
    part_cols += [_measurement_columns(pq.read_schema(path)) for path in paths[1:]]
    columns = [c for cols in part_cols for c in cols]
    values = np.empty((len(rows), len(columns)), dtype=np.float64)
    part_start = 0
    for index, (path, cols) in enumerate(zip(paths, part_cols)):
        handle = first_part if index == 0 else pq.ParquetFile(path)
        for start in range(0, len(cols), COLUMNS_PER_PART):
            block = handle.read(columns=cols[start : start + COLUMNS_PER_PART])
            for offset, column in enumerate(block.columns):
                values[:, part_start + start + offset] = column.to_numpy(zero_copy_only=False)[rows]
        part_start += len(cols)
    values[np.isnan(values)] = 0.0
    return meta[model_col].to_numpy()[rows], columns, values


def column_means(csv_path, columns):
    """The mean over all rows of each named measurement column, ignoring missing values, in the
    order given. Only the means are kept, so the table costs no more than a batch of it."""
    position = {c: k for k, c in enumerate(columns)}
    means = np.empty(len(columns))
    for path in _files(csv_path):
        wanted = [c for c in _measurement_columns(pq.read_schema(path)) if c in position]
        if not wanted:
            continue
        handle = pq.ParquetFile(path)
        for start in range(0, len(wanted), COLUMNS_PER_PART):
            batch = wanted[start : start + COLUMNS_PER_PART]
            for name, column in zip(batch, handle.read(columns=batch).columns):
                means[position[name]] = np.nanmean(column.to_numpy(zero_copy_only=False))
    return means
