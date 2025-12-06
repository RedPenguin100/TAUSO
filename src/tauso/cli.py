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
import time
from pyfaidx import Fasta
from gffutils.iterators import DataIterator
from tauso.data import get_paths, get_data_dir
from tauso.off_target.search import find_all_gene_off_targets, get_bowtie_index_base


def download_and_gunzip(url, dest_path):
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


# --- MAIN COMMAND ---

@click.group()
def main():
    """Tauso: ASO Design Toolkit"""
    pass


GENCODE_HUMAN_RELEASE = "34"


def get_genome_metadata(genome):
    """
    Returns download URLs for supported genomes.
    """
    CONFIG = {
        # --- MAMMALS (GENCODE) ---
        'GRCh38': {
            'base_url': f'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{GENCODE_HUMAN_RELEASE}',
            'fasta_name': 'GRCh38.primary_assembly.genome.fa.gz',
            'gtf_name': f'gencode.v{GENCODE_HUMAN_RELEASE}.basic.annotation.gtf.gz'
        },
        'GRCm39': {
            'base_url': 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33',
            'fasta_name': 'GRCm39.primary_assembly.genome.fa.gz',
            'gtf_name': 'gencode.vM33.basic.annotation.gtf.gz'
        },

        # --- YEAST (Ensembl) ---
        'R64-1-1': {
            'base_url': 'http://ftp.ensembl.org/pub/release-110',  # Stable release
            'fasta_name': 'fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz',
            'gtf_name': 'gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.110.gtf.gz'
        },

        # --- E. COLI (Ensembl Bacteria) ---
        # Note: E. coli K-12 MG1655
        'ASM584v2': {
            'base_url': 'http://ftp.ensemblgenomes.org/pub/bacteria/release-57',
            'fasta_name': 'fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.toplevel.fa.gz',
            'gtf_name': 'gtf/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.57.gtf.gz'
        }
    }

    if genome not in CONFIG:
        return None

    meta = CONFIG[genome]
    base = meta['base_url']

    # Logic to handle different URL structures between GENCODE and Ensembl
    # Ensembl URLs in config include the full relative path, GENCODE were simple appends
    # We can standardize by just joining base + name
    fasta_url = f"{base}/{meta['fasta_name']}"
    gtf_url = f"{base}/{meta['gtf_name']}"

    return fasta_url, gtf_url


@main.command()
@click.option('--genome', default='GRCh38', help='Genome name (GRCh38 or GRCm39).')
@click.option('--force', is_flag=True, help="Force re-download and rebuild.")
def setup_genome(genome, force):
    """
    Sets up the genome environment (Download -> Index -> Database).
    Supports GRCh38 (Human) and GRCm39 (Mouse).
    Always downloads the 'basic' gene annotation subset.

    :Note: You may set the TAUSO_DATA_DIR environment variable to switch the folder.
    """
    paths = get_paths(genome)
    fasta_path = paths['fasta']
    gtf_path = paths['gtf']
    db_path = paths['db']
    sentinel_path = db_path + ".success"

    if force:
        click.echo("Force flag detected. Cleaning up old files...")
        for f in [fasta_path, gtf_path, db_path, fasta_path + ".fai", sentinel_path]:
            if os.path.exists(f): os.remove(f)

    # --- PHASE 1: Download & Index ---
    try:
        # Check if we have automatic download support for this genome
        urls = get_genome_metadata(genome)

        if urls:
            fasta_url, gtf_url = urls

            # Only download if missing
            if not os.path.exists(fasta_path):
                click.echo(f"Downloading {genome} FASTA from GENCODE...")
                download_and_gunzip(fasta_url, fasta_path)

            if not os.path.exists(gtf_path):
                click.echo(f"Downloading {genome} GTF (Basic Annotation) from GENCODE...")
                download_and_gunzip(gtf_url, gtf_path)
        else:
            # Fallback for unsupported genomes
            if not os.path.exists(fasta_path) or not os.path.exists(gtf_path):
                click.echo(click.style(f"Error: Automatic download not supported for '{genome}'.", fg="red"))
                click.echo(f"Supported genomes: GRCh38, GRCm39")
                click.echo(f"For others, place {genome}.fa and {genome}.gtf manually in {get_data_dir()}")
                sys.exit(1)

        if not os.path.exists(fasta_path + ".fai"):
            click.echo("  Indexing FASTA file...")
            Fasta(fasta_path)
        click.echo("✓ Files ready.")

    except Exception as e:
        click.echo(click.style(f"Error during setup: {e}", fg="red"))
        sys.exit(1)

    # --- PHASE 3: Build Database ---
    # (Database building logic remains identical)
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

        # ... (Rest of existing DB build logic: introns, bulk loading, etc.) ...

        # Just to complete the block visually for you:
        print()  # Newline after progress bar

        # Calculate Unique Introns
        click.echo("  - Calculating unique introns...")
        intron_gen = db.create_introns()
        unique_introns = []
        seen_keys = set()
        for i, intron in enumerate(intron_gen):
            key = (intron.seqid, intron.start, intron.end, intron.strand)
            if key not in seen_keys:
                seen_keys.add(key)
                unique_introns.append(intron)
        del seen_keys

        # Bulk Load Introns
        click.echo(f"  - Bulk loading {len(unique_introns):,} introns...")
        del db
        db = gffutils.FeatureDB(db_path)
        # Set PRAGMA opts...
        with click.progressbar(length=len(unique_introns), label="    Importing") as bar:
            def monitor_gen():
                for x in unique_introns:
                    yield x
                    bar.update(1)

            db.update(monitor_gen(), merge_strategy='create_unique', disable_infer_genes=True,
                      disable_infer_transcripts=True)

        with open(sentinel_path, 'w') as f:
            f.write("Setup completed successfully.")

        click.echo(click.style(f"✓ Setup complete for {genome}.", fg="green"))

    except Exception as e:
        click.echo(click.style(f"\nError building database: {e}", fg="red"))
        if os.path.exists(db_path): os.remove(db_path)
        if os.path.exists(sentinel_path): os.remove(sentinel_path)
        sys.exit(1)

# --- NEW COMMAND: OFF-TARGET SEARCH ---
@main.command()
@click.argument('sequence')
@click.option('--genome', default='GRCh38', help='Genome name (default: GRCh38).')
@click.option('--mismatches', '-m', default=3, help='Max mismatches allowed.')
@click.option('--output', '-o', default=None, help='Output CSV file.')
def run_off_target(sequence, genome, mismatches, output):
    """
    Search for off-target binding sites for a given ASO sequence.
    """
    try:
        paths = get_paths(genome)
        # Construct the expected path to the "SUCCESS" sentinel file
        # This mirrors the logic in get_bowtie_index_base
        fasta_dir = os.path.dirname(paths['fasta'])
        sentinel_path = os.path.join(fasta_dir, f"{genome}_bowtie_index", "SUCCESS")

        if not os.path.exists(sentinel_path):
            click.echo(click.style(f"❌ Error: Bowtie index for '{genome}' is missing.", fg="red"))
            click.echo("The search cannot run because the genome index has not been built.")
            click.echo(f"\nPlease run this command first:\n  tauso setup-bowtie --genome {genome}")
            sys.exit(1)

    except Exception as e:
        # Catch issues like get_paths failing (e.g. if setup-genome wasn't run)
        click.echo(click.style(f"Error checking environment: {e}", fg="red"))
        sys.exit(1)

    original_seq = sequence
    sequence = sequence.upper().replace("U", "T")

    if original_seq != sequence:
        click.echo(f"Normalized input sequence: {original_seq} -> {sequence}")

    click.echo(f"Searching off-targets for {sequence} in {genome} (Max MM: {mismatches})...")

    start_time = time.time()
    try:
        hits_df = find_all_gene_off_targets(sequence, genome=genome, max_mismatches=mismatches)
    except Exception as e:
        click.echo(click.style(f"Search failed: {e}", fg="red"))
        sys.exit(1)

    duration = time.time() - start_time
    click.echo(f"Search completed in {duration:.2f}s. Found {len(hits_df)} hits.")

    if hits_df.empty:
        click.echo("No off-targets found.")
        return

    # Sort logic
    hits_df = hits_df.sort_values(by=['mismatches', 'chrom', 'start'])

    # Display Top 20
    click.echo("\nTop Hits:")
    click.echo(f"{'Chrom':<10} | {'Start':<12} | {'Str':<3} | {'MM':<2} | {'Gene':<15} | {'Region'}")
    click.echo("-" * 70)

    count = 0
    for _, hit in hits_df.iterrows():
        if count >= 20:
            click.echo(f"... and {len(hits_df) - 20} more.")
            break

        gene = str(hit['gene_name']) if hit['gene_name'] else "None"
        click.echo(
            f"{hit['chrom']:<10} | {hit['start']:<12} | {hit['strand']:<3} | {hit['mismatches']:<2} | {gene:<15} | {hit['region_type']}")
        count += 1

    if output:
        hits_df.to_csv(output, index=False)
        click.echo(f"\nFull results saved to: {output}")


@main.command()
@click.option('--genome', default='GRCh38', help='Genome name (default: GRCh38).')
@click.option('--force', is_flag=True, help="Force rebuild of the index.")
def setup_bowtie(genome, force):
    """
    Generates (or validates) the Bowtie index for the specified genome.
    Requires 'setup-genome' to be run first to ensure the FASTA file exists.
    """
    click.echo(f"Initializing Bowtie setup for {genome}...")

    try:
        # We assume the FASTA exists. If not, get_bowtie_index_base throws FileNotFoundError
        # This function handles the logic:
        # 1. Check if index exists -> return if so.
        # 2. If --force or missing -> run bowtie-build.
        index_path = get_bowtie_index_base(genome=genome, force_rebuild=force)

        click.echo(click.style(f"✓ Bowtie index ready at: {index_path}", fg="green"))

    except FileNotFoundError:
        click.echo(click.style(f"Error: FASTA file for {genome} not found.", fg="red"))
        click.echo(f"Please run 'tauso setup-genome --genome {genome}' first.")
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"Error building Bowtie index: {e}", fg="red"))
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