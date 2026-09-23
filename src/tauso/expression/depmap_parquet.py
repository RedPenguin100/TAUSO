"""DepMap expression tables as Parquet, converted and read without knowing how they are stored.

Every function takes the table's CSV path (the name DepMap publishes it under), whether or not
the CSV is still there. The CSVs are tens of thousands of columns wide, and a Parquet writer's
memory follows the width it is given, so a table is stored as a directory of narrow Parquet
parts: the first holds the profile columns, each later one COLUMNS_PER_PART expression columns,
and every part the same rows in the same order. A data dir converted before the split holds one
Parquet file instead, which reads the same way.
"""

import hashlib
import os
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from ..cli_utils import sha256_file

DEPMAP_PROFILE_COLUMNS = frozenset(
    {"Unnamed: 0", "", "SequencingID", "ModelID", "IsDefaultEntryForModel", "ModelConditionID", "IsDefaultEntryForMC"}
)
"""The columns of a DepMap expression table that describe the sequencing run, not a measurement."""

COLUMNS_PER_PART = 2000
"""Expression columns per Parquet part, and per read. A Parquet writer holds buffers for every
column until a row group is written, so narrow parts keep its memory small."""


def _parts_dir(csv_path):
    path = Path(csv_path)
    return path.with_name(path.stem + "_parts")


def _single_file(csv_path):
    return Path(csv_path).with_suffix(".parquet")


def _location(csv_path):
    """Where the Parquet is: the parts directory, or the single file of an older data dir."""
    single = _single_file(csv_path)
    return single if single.exists() and not _parts_dir(csv_path).is_dir() else _parts_dir(csv_path)


def _files(csv_path):
    """The Parquet files in column order; empty if the CSV has not been converted."""
    location = _location(csv_path)
    if location.is_dir():
        return sorted(str(p) for p in location.glob("part-*.parquet"))
    return [str(location)] if location.exists() else []


def _sha256_file_path(csv_path):
    location = _location(csv_path)
    return location.with_name(location.name + ".sha256")


def _expression_columns_in(schema):
    """All but the profile columns and a stored pandas index."""
    index_cols = {c for c in (schema.pandas_metadata or {}).get("index_columns", []) if isinstance(c, str)}
    return [c for c in schema.names if c not in DEPMAP_PROFILE_COLUMNS and c not in index_cols]


def parquet_exists(csv_path):
    """Whether the CSV has been converted to Parquet."""
    return bool(_files(csv_path))


def expression_columns(csv_path):
    """The table's gene or transcript columns, in order."""
    return [c for path in _files(csv_path) for c in _expression_columns_in(pq.read_schema(path))]


def _split_csv_by_columns(csv_path, out_dir):
    """Split the CSV into narrow CSVs of the same rows, and return their paths in column order.

    The split is plain text, one line at a time, so no parser ever sees the full width.
    """
    with open(csv_path) as source:
        header = source.readline().rstrip("\r\n").split(",")
        profile = [i for i, c in enumerate(header) if c in DEPMAP_PROFILE_COLUMNS]
        expression = [i for i, c in enumerate(header) if c not in DEPMAP_PROFILE_COLUMNS]
        groups = [profile] + [
            expression[start : start + COLUMNS_PER_PART] for start in range(0, len(expression), COLUMNS_PER_PART)
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


def remove_parquet(csv_path):
    """Remove the CSV's Parquet and its saved hash, however they are stored."""
    for location in (_parts_dir(csv_path), _single_file(csv_path)):
        if location.is_dir():
            shutil.rmtree(location)
        else:
            location.unlink(missing_ok=True)
        location.with_name(location.name + ".sha256").unlink(missing_ok=True)


def csv_to_parquet(csv_path):
    """Convert the CSV to Parquet, replacing any earlier conversion, save its hash, and remove the CSV.

    Each narrow part is parsed by pandas, which measured lighter than pyarrow on both tables.
    """
    directory = _parts_dir(csv_path)
    partial = directory.with_name(directory.name + ".partial")
    shutil.rmtree(partial, ignore_errors=True)
    partial.mkdir()
    for part_csv in _split_csv_by_columns(csv_path, partial):
        pd.read_csv(part_csv).to_parquet(part_csv.with_suffix(".parquet"), index=False)
        part_csv.unlink()
    remove_parquet(csv_path)
    # Named only once it is whole, so a killed conversion is not mistaken for a finished one.
    os.replace(partial, directory)
    os.remove(csv_path)
    save_parquet_sha256(csv_path)


def parquet_sha256(csv_path):
    """The Parquet's SHA-256: of its single file, or of its parts' names and hashes."""
    location = _location(csv_path)
    if not location.is_dir():
        return sha256_file(str(location))
    listing = "".join(f"{Path(p).name} {sha256_file(p)}\n" for p in _files(csv_path))
    return hashlib.sha256(listing.encode()).hexdigest()


def saved_parquet_sha256(csv_path):
    """The hash saved when the Parquet was written, or None if none was."""
    path = _sha256_file_path(csv_path)
    return path.read_text().strip() if path.exists() else None


def save_parquet_sha256(csv_path):
    """Save the Parquet's hash beside it, so a later run can tell it has changed on disk."""
    _sha256_file_path(csv_path).write_text(parquet_sha256(csv_path))


def _default_rows(first_part, model_ids):
    """The rows to read for the wanted models, in table order, and the model id of each.

    A model can carry several sequencing profiles; the one DepMap flags as the default is read,
    and a model with two left is read from its first.
    """
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
    return rows, meta[model_col].to_numpy()[rows]


def _column_blocks(paths, first_part, part_cols, rows):
    """Yield (index of its first column, the rows' values) for each block of up to COLUMNS_PER_PART
    columns, in column order, missing values as 0. One part is open at a time: every open part
    holds its footer in memory."""
    part_start = 0
    for index, (path, cols) in enumerate(zip(paths, part_cols)):
        handle = first_part if index == 0 else pq.ParquetFile(path)
        for start in range(0, len(cols), COLUMNS_PER_PART):
            batch = cols[start : start + COLUMNS_PER_PART]
            block = np.empty((len(rows), len(batch)), dtype=np.float64)
            for offset, column in enumerate(handle.read(columns=batch).columns):
                block[:, offset] = column.to_numpy(zero_copy_only=False)[rows]
            block[np.isnan(block)] = 0.0
            yield part_start + start, block
        part_start += len(cols)


def _wanted_blocks(csv_path, model_ids):
    """The wanted cell lines' model ids, the expression columns, and their values as column blocks.

    Parquet stores a table by column, so no cell line can be read without decoding every column
    whole: the blocks come one part at a time, each holding every wanted cell line.
    """
    paths = _files(csv_path)
    first_part = pq.ParquetFile(paths[0])
    rows, found = _default_rows(first_part, model_ids)
    part_cols = [_expression_columns_in(first_part.schema_arrow)]
    part_cols += [_expression_columns_in(pq.read_schema(path)) for path in paths[1:]]
    columns = [c for cols in part_cols for c in cols]
    return found, columns, _column_blocks(paths, first_part, part_cols, rows)


def read_cell_lines(csv_path, model_ids):
    """The wanted cell lines' expression, as (model ids found, expression columns, values).

    The default profile of each, missing values as 0, all held at once: for a cohort small
    enough to hold. iter_cell_lines gives the same one cell line at a time.
    """
    found, columns, blocks = _wanted_blocks(csv_path, model_ids)
    values = np.empty((len(found), len(columns)), dtype=np.float64)
    for first_column, block in blocks:
        values[:, first_column : first_column + block.shape[1]] = block
    return found, columns, values


def iter_cell_lines(csv_path, model_ids):
    """Yield (model id, expression columns, values) for each wanted cell line, as read_cell_lines
    reads them, in table order.

    A cell line is only whole once every part has been decoded, so the blocks are written to a
    temporary file, one cell line per row, and the cell lines are read back from it one at a time.
    What stays in memory is one block and one cell line however large the cohort; the file, in
    the data dir because some machines keep the system temp dir in memory, is 8 bytes per value
    (1.9 MB per cell line of the transcript table) and is gone when the last one is yielded.
    """
    found, columns, blocks = _wanted_blocks(csv_path, model_ids)
    width = len(columns)
    with tempfile.TemporaryFile(dir=Path(csv_path).parent) as scratch:
        for first_column, block in blocks:
            for position, values in enumerate(block):
                scratch.seek((position * width + first_column) * 8)
                scratch.write(values.tobytes())
        for position, model_id in enumerate(found):
            scratch.seek(position * width * 8)
            yield model_id, columns, np.fromfile(scratch, dtype=np.float64, count=width)


def mean_expression(csv_path, columns):
    """Each named column's mean over every cell line, ignoring missing values, in the order given.

    Only the means are kept, so the table costs no more than a batch of it.
    """
    position = {c: k for k, c in enumerate(columns)}
    means = np.empty(len(columns))
    for path in _files(csv_path):
        wanted = [c for c in _expression_columns_in(pq.read_schema(path)) if c in position]
        if not wanted:
            continue
        handle = pq.ParquetFile(path)
        for start in range(0, len(wanted), COLUMNS_PER_PART):
            batch = wanted[start : start + COLUMNS_PER_PART]
            for name, column in zip(batch, handle.read(columns=batch).columns):
                means[position[name]] = np.nanmean(column.to_numpy(zero_copy_only=False))
    return means
