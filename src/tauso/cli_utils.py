"""Shared helpers for tauso CLI commands."""

import gzip
import hashlib
import itertools
import os
import shutil
import sys

import click
import requests


def sha1_file(path: str) -> str:
    h = hashlib.sha1()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def md5_file(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def echo_err(msg: str) -> None:
    click.echo(click.style(f"❌ {msg}", fg="red"))


def echo_warn(msg: str) -> None:
    click.echo(click.style(f"⚠ {msg}", fg="yellow"))


def echo_ok(msg: str) -> None:
    click.echo(click.style(f"✓ {msg}", fg="green"))


def verify_hash_or_exit(path: str, expected: str, algo: str = "sha1") -> None:
    """Compute a hash of `path` and exit nonzero if it does not match `expected`."""
    hashers = {"sha1": sha1_file, "sha256": sha256_file, "md5": md5_file}
    if algo not in hashers:
        raise ValueError(f"Unsupported hash algorithm: {algo}")
    actual = hashers[algo](path)
    if actual != expected:
        echo_err(f"{algo.upper()} mismatch for {os.path.basename(path)}")
        click.echo(f"  Expected: {expected}")
        click.echo(f"  Got:      {actual}")
        sys.exit(1)


def download_with_progress(url: str, dest: str, label: str | None = None) -> None:
    """Stream `url` to `dest` with a click progress bar. Cleans up `dest` on failure."""
    label = label or f"Downloading {os.path.basename(dest)}"
    try:
        with requests.get(url, stream=True, timeout=30) as r:
            r.raise_for_status()
            total = int(r.headers.get("content-length", 0))
            with open(dest, "wb") as f, click.progressbar(length=total, label=label) as bar:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
                    bar.update(len(chunk))
    except Exception:
        if os.path.exists(dest):
            os.remove(dest)
        raise


def download_and_gunzip(url: str, dest_path: str, remove_gz: bool = False) -> None:
    """Download `url` (gzipped) and decompress to `dest_path`."""
    if os.path.exists(dest_path):
        click.echo(f"  File already exists: {os.path.basename(dest_path)}")
        return

    temp_gz = dest_path + ".gz"
    try:
        download_with_progress(url, temp_gz, label=f"    Downloading {os.path.basename(url)}")

        click.echo(f"  Unzipping to {os.path.basename(dest_path)}...")
        with gzip.open(temp_gz, "rb") as f_in, open(dest_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        if remove_gz and os.path.exists(temp_gz):
            os.remove(temp_gz)
    except Exception:
        for p in (dest_path, temp_gz):
            if os.path.exists(p):
                os.remove(p)
        raise


def count_lines(filepath: str) -> int:
    with open(filepath, "rb") as f:
        return sum(1 for line in f if not line.startswith(b"#"))


def batch_iterator(iterator, batch_size: int = 1000):
    while True:
        batch = list(itertools.islice(iterator, batch_size))
        if not batch:
            break
        yield batch
