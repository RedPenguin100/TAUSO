import os

import click

from tauso.cli_utils import download_with_progress, echo_ok, echo_warn, sha1_file, verify_hash_or_exit

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
