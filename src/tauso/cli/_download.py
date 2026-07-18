import os
import sys

import click
import pandas as pd

from tauso.cli_utils import download_with_progress, echo_err, echo_ok, echo_warn, sha1_file, verify_hash_or_exit

# Zenodo-mirrored DepMap "Public 25Q3" snapshot. DepMap silently re-uploads files
# under the same release name, so we pin to an immutable Zenodo record and verify
# the SHA1 of each download.
ZENODO_DEPMAP_RECORD = "20355477"
DEPMAP_FILES_SHA1 = {
    "Model.csv": "4e9805ecf79d187e1fb5d4c760312e5a40729e34",
    "OmicsProfiles.csv": "fc5a1ed86ea89f805d56715f439e9738b3e28a72",
    "OmicsExpressionTPMLogp1HumanAllGenesStranded.csv": "22ac03aa45a6b9ef4f60e9ed8bb574e64dcb56f6",
}

# rRNA off-target reference: human cytoplasmic 18S/5.8S/28S/5S RefSeq records, frozen on Zenodo.
ZENODO_RRNA_RECORD = "21071791"
RRNA_SHA1 = "377ded75d51e57a12eee899c333f30431de0f7ff"


def _zenodo_file_url(record_id: str, filename: str) -> str:
    return f"https://zenodo.org/records/{record_id}/files/{filename}"


def _zenodo_content_url(record_id: str, filename: str) -> str:
    return f"https://zenodo.org/api/records/{record_id}/files/{filename}/content"


def _ensure_depmap_file(filename: str, expected_sha1: str, data_dir: str, force: bool) -> bool:
    """Ensure `filename` exists in `data_dir` with the pinned SHA1. Returns True if the file
    was (re-)downloaded, False if an existing valid copy was reused."""
    dest = os.path.join(data_dir, filename)

    if os.path.exists(dest) and not force:
        if sha1_file(dest) == expected_sha1:
            echo_ok(f"{filename} exists (SHA1 verified).")
            return False
        echo_warn(f"SHA1 mismatch for {filename} — re-downloading.")

    url = _zenodo_file_url(ZENODO_DEPMAP_RECORD, filename)
    click.echo(f"Downloading {filename} from Zenodo...")
    download_with_progress(url, dest, label=f"    {filename}")
    verify_hash_or_exit(dest, expected_sha1, algo="sha1")
    echo_ok(f"Downloaded {filename} (SHA1 verified).")
    return True


def _ensure_zenodo_content_file(
    record_id: str, filename: str, destination: str, expected_hash: str, algo: str, force: bool
) -> None:
    """Ensure `destination` holds Zenodo file `filename`, verified against `expected_hash`."""
    if os.path.exists(destination) and not force:
        verify_hash_or_exit(destination, expected_hash, algo=algo)
        echo_ok(f"Existing {filename} matches expected {algo.upper()}. Skipping download.")
        return
    try:
        download_with_progress(_zenodo_content_url(record_id, filename), destination, label=f"Downloading {filename}")
        verify_hash_or_exit(destination, expected_hash, algo=algo)
        echo_ok(f"Downloaded and verified: {destination}")
    except Exception as e:
        echo_err(f"Error downloading {filename}: {e}")
        sys.exit(1)


def _ensure_zenodo_table(
    record_id: str,
    source_name: str,
    source_path: str,
    parquet_path: str,
    expected_hash: str,
    algo: str,
    force: bool,
    usecols=None,
) -> bool:
    """Download a tabular Zenodo file to `source_path`, verify its hash, convert it to
    `parquet_path` (optionally keeping only `usecols`), and drop the source. Returns False
    (skips) if the parquet already exists and not `force`."""
    if os.path.exists(parquet_path) and not force:
        echo_ok(f"{os.path.basename(parquet_path)} already present. Skipping.")
        return False
    try:
        download_with_progress(
            _zenodo_content_url(record_id, source_name), source_path, label=f"Downloading {source_name}"
        )
        verify_hash_or_exit(source_path, expected_hash, algo=algo)
        echo_ok(f"Downloaded and verified: {source_path}")
        click.echo("  Converting to Parquet...")
        pd.read_csv(source_path, usecols=usecols).to_parquet(parquet_path, index=False)
        echo_ok(f"Converted to Parquet: {parquet_path}")
        os.remove(source_path)
        return True
    except Exception as e:
        echo_err(f"Error setting up {os.path.basename(parquet_path)}: {e}")
        sys.exit(1)
