import os
import subprocess
import sys
import gzip
import shutil
from importlib.resources import files

import requests
import click
import gffutils
import itertools
from platformdirs import user_data_dir
from pyfaidx import Fasta
from gffutils.iterators import DataIterator

# --- CONFIGURATION ---
APP_NAME = "tauso"
GENCODE_RELEASE = "49"
GENCODE_BASE_URL = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{GENCODE_RELEASE}"


def get_data_dir():
    """Returns the path to ~/.local/share/tauso/ (or OS equivalent)."""
    d = user_data_dir(APP_NAME)
    os.makedirs(d, exist_ok=True)
    return d


def get_gencode_urls(subset):
    """Returns URLs for the PRIMARY ASSEMBLY."""
    fasta_url = f"{GENCODE_BASE_URL}/GRCh38.primary_assembly.genome.fa.gz"
    if subset == 'basic':
        gtf_name = f"gencode.v{GENCODE_RELEASE}.primary_assembly.basic.annotation.gtf.gz"
    else:
        gtf_name = f"gencode.v{GENCODE_RELEASE}.primary_assembly.annotation.gtf.gz"
    gtf_url = f"{GENCODE_BASE_URL}/{gtf_name}"
    return fasta_url, gtf_url


def download_and_gunzip(url, dest_path):
    """Downloads and unzips a file with a progress bar."""
    if os.path.exists(dest_path):
        click.echo(f"  File already exists: {os.path.basename(dest_path)}")
        return

    try:
        click.echo(f"  Downloading {os.path.basename(url)}...")
        temp_gz = dest_path + ".gz"

        with requests.get(url, stream=True, timeout=30) as r:
            r.raise_for_status()
            total_size = int(r.headers.get('content-length', 0))

            with open(temp_gz, 'wb') as f:
                with click.progressbar(length=total_size, label="    Downloading") as bar:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
                        bar.update(len(chunk))

        click.echo(f"  Unzipping to {os.path.basename(dest_path)}...")
        with gzip.open(temp_gz, 'rb') as f_in:
            with open(dest_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        os.remove(temp_gz)
    except Exception as e:
        if os.path.exists(dest_path): os.remove(dest_path)
        if os.path.exists(temp_gz): os.remove(temp_gz)
        raise e


def count_lines(filepath):
    with open(filepath, 'rb') as f:
        return sum(1 for _ in f)


def batch_iterator(iterator, batch_size=1000):
    while True:
        batch = list(itertools.islice(iterator, batch_size))
        if not batch: break
        yield batch


# --- CLI COMMAND ---

@click.group()
def main():
    pass


@main.command()
@click.option('--subset', type=click.Choice(['basic', 'full']), default='basic',
              help='Annotation subset: "basic" (representative) or "full" (comprehensive).')
@click.option('--force', is_flag=True, help="Force re-download and rebuild.")
def setup_genome(subset, force):
    """
    Sets up the genome environment (Download -> Index -> Database).
    """
    data_dir = get_data_dir()
    fasta_path = os.path.join(data_dir, "GRCh38.fa")
    gtf_path = os.path.join(data_dir, "GRCh38.gtf")
    db_path = os.path.join(data_dir, "GRCh38.db")
    sentinel_path = db_path + ".success"

    if force:
        click.echo("Force flag detected. Cleaning up old files...")
        for f in [fasta_path, gtf_path, db_path, fasta_path + ".fai", sentinel_path]:
            if os.path.exists(f): os.remove(f)

    # --- PHASE 1: Download Raw Files ---
    try:
        fasta_url, gtf_url = get_gencode_urls(subset)
        click.echo(f"Setting up in: {data_dir}")

        # Only download if missing
        download_and_gunzip(fasta_url, fasta_path)
        download_and_gunzip(gtf_url, gtf_path)

    except Exception as e:
        click.echo(click.style(f"Error during download: {e}", fg="red"))
        sys.exit(1)

    # --- PHASE 2: Index FASTA ---
    # This creates the .fai file allowing instant random access to sequences
    try:
        if not os.path.exists(fasta_path + ".fai"):
            click.echo("  Indexing FASTA file (this happens once)...")
            Fasta(fasta_path)
        click.echo("✓ Files ready.")

    except Exception as e:
        click.echo(click.style(f"Error during indexing: {e}", fg="red"))
        sys.exit(1)

    # --- PHASE 3: Build Database ---
    if os.path.exists(db_path):
        if os.path.exists(sentinel_path):
            click.echo(f"✓ Database already exists at {db_path}")
            return
        else:
            click.echo(click.style(f"⚠ Found incomplete/corrupt database. Rebuilding...", fg="yellow"))
            if os.path.exists(db_path): os.remove(db_path)

    try:
        click.echo(f"Building database at {db_path}...")

        # 1. Parsing GTF
        total_lines = count_lines(gtf_path)
        data_it = DataIterator(gtf_path)

        click.echo(f"  - Parsing {total_lines:,} annotation lines...")
        with click.progressbar(length=total_lines, label="    Parsing GTF") as bar:
            def progress_wrapper():
                for feature in data_it:
                    bar.update(1)
                    yield feature

            db = gffutils.create_db(
                progress_wrapper(),
                dbfn=db_path,
                force=True,
                keep_order=True,
                merge_strategy='merge',
                sort_attribute_values=True,
                disable_infer_genes=True,
                disable_infer_transcripts=True
            )

        # 2. Calculating & Deduplicating Introns (IN RAM)
        click.echo("  - Calculating unique introns...")
        intron_gen = db.create_introns()

        unique_introns = []
        seen_keys = set()

        for i, intron in enumerate(intron_gen):
            key = (intron.seqid, intron.start, intron.end, intron.strand)
            if key not in seen_keys:
                seen_keys.add(key)
                unique_introns.append(intron)

            if i % 50000 == 0:
                sys.stdout.write(f"\r    Scanned {i:,} isoforms -> Kept {len(unique_introns):,} unique...")
                sys.stdout.flush()

        click.echo(f"\r    Scanned {i:,} isoforms -> Kept {len(unique_introns):,} unique. Done.")
        del seen_keys

        # 3. Bulk Loading (RAM -> DB)
        click.echo(f"  - Bulk loading {len(unique_introns):,} introns...")

        del db
        db = gffutils.FeatureDB(db_path)

        # MAX SPEED SETTINGS
        db.execute("PRAGMA synchronous = OFF")
        db.execute("PRAGMA journal_mode = MEMORY")
        db.execute("PRAGMA cache_size = 1000000")
        db.execute("PRAGMA locking_mode = EXCLUSIVE")

        with click.progressbar(length=len(unique_introns), label="    Importing") as bar:
            def monitor_gen():
                for x in unique_introns:
                    yield x
                    bar.update(1)

            db.update(monitor_gen(),
                      merge_strategy='create_unique',
                      disable_infer_genes=True,
                      disable_infer_transcripts=True)

        print()

        # Mark Success
        with open(sentinel_path, 'w') as f:
            f.write("Setup completed successfully.")

        click.echo(click.style(f"✓ Setup complete.", fg="green"))

    except Exception as e:
        click.echo(click.style(f"\nError building database: {e}", fg="red"))
        if os.path.exists(db_path): os.remove(db_path)
        if os.path.exists(sentinel_path): os.remove(sentinel_path)
        sys.exit(1)

@main.command(context_settings=dict(
    ignore_unknown_options=True,
    allow_extra_args=True
))
@click.pass_context
@click.option(
    "--force-clone", "-f",
    is_flag=True,
    help="Force re-cloning of the raccess repository (passed to install_raccess.sh).",
)
def install_raccess(ctx, force_clone):
    """
    Run the raccess installation per their license
    """
    script_path = files('tauso._raccess') / 'install_raccess.sh'

    forwarded_args = list(ctx.args)
    if force_clone:
        forwarded_args.append("--force-clone")  # or "-f", your choice

    cmd = ['bash', str(script_path)] + forwarded_args

    print(f"+ {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


if __name__ == '__main__':
    main()
