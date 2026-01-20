import zipfile
from pathlib import Path

import pandas as pd
import json
import re
import numpy as np

import io
import pandas as pd
import requests
import click
import os
import subprocess
import sys
import gdown
import gzip
import shutil
from importlib.resources import files

import gffutils
import itertools
import time
from pyfaidx import Fasta
from gffutils.iterators import DataIterator

from tauso.features.codon_usage.find_cai_reference import load_cell_line_gene_maps
from tauso.data.data import get_paths, get_data_dir
from tauso.features.cai import calc_CAI_weight
from tauso.genome.TranscriptMapper import GeneCoordinateMapper, build_gene_sequence_registry
from tauso.genome.read_human_genome import get_locus_to_data_dict
from tauso.off_target.search import find_all_gene_off_targets, get_bowtie_index_base

import click


# --- MAIN COMMAND ---

@click.group()
def main():
    """Tauso: ASO Design Toolkit"""
    pass


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


@main.command()
@click.option('--force', is_flag=True, help="Force redownload.")
def setup_attract(force):
    """
    Downloads ATtRACT database.
    1. Extracts 'pwm.txt' (The matrices).
    2. Filters 'ATtRACT_db.txt' for Homo sapiens and saves as 'RBS_motifs_Homo_sapiens.csv'.
    """
    data_dir = get_data_dir()
    os.makedirs(data_dir, exist_ok=True)

    # Define output paths
    metadata_output = os.path.join(data_dir, "RBS_motifs_Homo_sapiens.csv")
    pwm_output = os.path.join(data_dir, "pwm.txt")

    if os.path.exists(metadata_output) and os.path.exists(pwm_output) and not force:
        click.echo(f"✓ ATtRACT data already exists in {data_dir}")
        return

    url = "https://attract.cnic.es/attract/static/ATtRACT.zip"
    click.echo("Downloading ATtRACT database...")

    try:
        r = requests.get(url)
        r.raise_for_status()

        with zipfile.ZipFile(io.BytesIO(r.content)) as z:
            # 1. Find and Extract pwm.txt
            try:
                pwm_filename = next(name for name in z.namelist() if name.endswith('pwm.txt'))
                click.echo(f"  Extracting {pwm_filename}...")
                with z.open(pwm_filename) as source, open(pwm_output, "wb") as target:
                    target.write(source.read())
            except StopIteration:
                raise FileNotFoundError("pwm.txt not found in the archive.")

            # 2. Find and Process ATtRACT_db.txt
            try:
                db_filename = next(name for name in z.namelist() if name.endswith('ATtRACT_db.txt'))
            except StopIteration:
                raise FileNotFoundError("ATtRACT_db.txt not found in the archive.")

            click.echo(f"  Processing {db_filename}...")
            with z.open(db_filename) as f:
                # Read original tab-separated file
                df = pd.read_csv(f, sep='\t')

                click.echo(f"  Filtering {len(df)} entries for Homo sapiens...")
                human_df = df[df['Organism'] == 'Homo_sapiens'].copy()

                if human_df.empty:
                    raise ValueError("No Homo_sapiens entries found.")

                # Save as standard CSV
                human_df.to_csv(metadata_output, index=False)
                click.echo(f"✓ Saved metadata to {metadata_output}")
                click.echo(f"✓ Saved matrices to {pwm_output}")

    except Exception as e:
        click.echo(click.style(f"❌ Error: {e}", fg="red"))
        sys.exit(1)


@main.command()
@click.option('--force', is_flag=True, help="Force redownload.")
def setup_depmap(force):
    """
    Downloads DepMap Public 25Q3 data directly from the DepMap API.
    Auto-fetches fresh signed URLs to ensure downloads work.
    """
    data_dir = get_data_dir()
    os.makedirs(data_dir, exist_ok=True)

    # 1. Configuration
    RELEASE = "DepMap Public 25Q3"
    TARGET_FILES = {
        # Remote Filename -> Local Filename (can be same)
        "Model.csv": "Model.csv",  # <--- Added Model.csv here
        "OmicsProfiles.csv": "OmicsProfiles.csv",
        "OmicsExpressionTPMLogp1HumanAllGenesStranded.csv": "OmicsExpressionTPMLogp1HumanAllGenesStranded.csv"
    }

    click.echo(f"Initializing DepMap setup for: {RELEASE}")
    click.echo("Fetching fresh download URLs from DepMap API...")

    # 2. Fetch the Master Index CSV
    # This endpoint returns a CSV with columns: release, filename, url, etc.
    index_url = "https://depmap.org/portal/api/download/files"

    try:
        r = requests.get(index_url)
        r.raise_for_status()
        # Parse directly into Pandas
        df = pd.read_csv(io.StringIO(r.text))
    except Exception as e:
        click.echo(click.style(f"❌ Error fetching DepMap index: {e}", fg="red"))
        sys.exit(1)

    # 3. Filter for 25Q3
    # Note: Column names are 'release', 'filename', 'url'
    release_df = df[df['release'] == RELEASE]

    if release_df.empty:
        click.echo(click.style(f"❌ Release '{RELEASE}' not found in API.", fg="red"))
        click.echo(f"Available releases: {df['release'].unique()[:5]}...")
        sys.exit(1)

    # 4. Download Loop
    for remote_name, local_name in TARGET_FILES.items():
        dest = os.path.join(data_dir, local_name)

        if os.path.exists(dest) and not force:
            click.echo(f"✓ {local_name} exists.")
            continue

        # Find the row for this specific file
        file_row = release_df[release_df['filename'] == remote_name]

        if file_row.empty:
            click.echo(click.style(f"⚠ Warning: File '{remote_name}' not found in {RELEASE}.", fg="yellow"))
            continue

        # Extract the signed URL
        download_url = file_row.iloc[0]['url']

        click.echo(f"Downloading {local_name}...")
        try:
            # Using wget with -O to save to destination
            subprocess.run(["wget", "-q", "--show-progress", "-O", dest, download_url], check=True)
            click.echo(click.style(f"✓ Downloaded {local_name}", fg="green"))
        except subprocess.CalledProcessError:
            click.echo(click.style(f"❌ Failed to download {local_name}", fg="red"))
            if os.path.exists(dest): os.remove(dest)

    click.echo("\nDepMap setup complete.")



@main.command()
@click.argument('cell_names', nargs=-1)
@click.option('--reset', is_flag=True, help="Clear existing list before adding.")
def add_cell(cell_names, reset):
    """
    Search for cell lines by name and add them to the project cohort.
    Example: tauso add-cell "HepG2" "HeLa" "U251"
    """
    data_dir = get_data_dir()
    profiles_path = os.path.join(data_dir, "OmicsProfiles.csv")
    manifest_path = os.path.join(data_dir, "cell_cohort.json")

    if not os.path.exists(profiles_path):
        click.echo(click.style("Error: OmicsProfiles.csv not found. Run 'setup-depmap' first.", fg="red"))
        sys.exit(1)

    # 1. Load Existing Manifest
    cohort = {}
    if os.path.exists(manifest_path) and not reset:
        with open(manifest_path, 'r') as f:
            cohort = json.load(f)

    # 2. Load Metadata (Optimized)
    click.echo("Loading metadata...")
    df = pd.read_csv(profiles_path, usecols=['ModelID', 'StrippedCellLineName', 'IsDefaultEntryForModel'])
    # Normalize for fuzzy search
    df['clean'] = df['StrippedCellLineName'].str.replace(r'[^a-zA-Z0-9]', '', regex=True).str.upper()

    # 3. Search
    for query in cell_names:
        clean_q = query.replace('-', '').replace(' ', '').upper()

        # Exact match
        match = df[df['clean'] == clean_q]
        # Fuzzy match
        if match.empty:
            match = df[df['clean'].str.contains(clean_q, na=False)]

        if not match.empty:
            # Prefer default entry
            best = match[match['IsDefaultEntryForModel'] == True]
            if best.empty: best = match.iloc[[0]]

            found_name = best.iloc[0]['StrippedCellLineName']
            ach_id = best.iloc[0]['ModelID']

            cohort[found_name] = ach_id
            click.echo(click.style(f"✓ Found: {query} -> {found_name} ({ach_id})", fg="green"))
        else:
            click.echo(click.style(f"⚠ Not Found: {query}", fg="yellow"))

    # 4. Save
    with open(manifest_path, 'w') as f:
        json.dump(cohort, f, indent=4)

    click.echo(f"Cohort saved to {manifest_path} ({len(cohort)} cell lines).")

@main.command()
@click.option('--genome', default='GRCh38', help='Genome version (default: GRCh38).')
def build_omics(genome):
    """
    Generates gene-level expression files for all cell lines in the cohort.
    Uses DepMap 'AllGenes' format (Gene Level Log2(TPM+1)).
    """
    if genome != 'GRCh38':
        # While the logic is genome-agnostic now, we keep this guard if your downstream tools expect GRCh38
        click.echo(click.style(f"⚠ Warning: This command is optimized for DepMap (Human GRCh38) data.", fg="yellow"))

    data_dir = get_data_dir()
    manifest_path = os.path.join(data_dir, "cell_cohort.json")
    exp_path = os.path.join(data_dir, "OmicsExpressionTPMLogp1HumanAllGenesStranded.csv")

    if not os.path.exists(manifest_path):
        click.echo(click.style("❌ No cohort found. Use 'tauso add-cell' first.", fg="red"))
        return

    if not os.path.exists(exp_path):
        click.echo(click.style(f"❌ Expression file not found: {exp_path}", fg="red"))
        click.echo("Run 'tauso setup-depmap' to download it.")
        return

    with open(manifest_path, 'r') as f:
        cohort = json.load(f)

    # We need the set of ACH-IDs to find in the huge CSV
    target_ids = set(cohort.values())
    click.echo(f"Processing {len(target_ids)} cell lines from cohort...")

    # --- Step 1: Parse Header & Clean Gene Names ---
    click.echo(f"Scanning {os.path.basename(exp_path)}...")

    with open(exp_path, 'r') as f:
        header_line = f.readline().strip()
        header = header_line.split(',')

    try:
        model_idx = header.index('ModelID')
    except ValueError:
        # DepMap CSVs usually have ModelID as the first column (index 0) if not named explicitly
        model_idx = 0

        # Identify Gene Columns: All columns AFTER the metadata (usually everything after ModelID)
    # The header looks like: ,ENSG000..., RPH3AL-AS1 (100506388), ...
    # We skip the first column (empty in your snippet) and ModelID column.

    # We'll assume any column that isn't ModelID is a data column.
    # Let's dynamically find the start of data.
    data_indices = [i for i, c in enumerate(header) if i != model_idx and c.strip()]

    # Pre-calculate clean gene names for the header
    # Format 1: "RPH3AL-AS1 (100506388)" -> "RPH3AL-AS1"
    # Format 2: "ENSG00000262038.1" -> "ENSG00000262038.1" (Keep as is, or strip version)
    clean_genes = []
    gene_regex = re.compile(r"^(.+?) \(\d+\)$")  # Capture "Symbol" from "Symbol (ID)"

    for i in data_indices:
        raw_col = header[i]
        match = gene_regex.match(raw_col)
        if match:
            clean_genes.append(match.group(1))
        else:
            clean_genes.append(raw_col)

    # --- Step 2: Stream & Extract ---
    output_dir = os.path.join(data_dir, "processed_expression")
    os.makedirs(output_dir, exist_ok=True)

    found_count = 0

    with open(exp_path, 'r') as f:
        # Skip header since we read it
        next(f)

        for line in f:
            # Optimization: Only split the start of the line to check ModelID
            # This avoids splitting 20k+ columns for rows we don't need
            row_start = line.split(',', model_idx + 1)
            if len(row_start) <= model_idx: continue

            curr_id = row_start[model_idx]

            if curr_id in target_ids:
                click.echo(f"  Extracting {curr_id}...")

                # Now split the full line
                parts = line.strip().split(',')

                # Extract values corresponding to our data indices
                vals = []
                for idx in data_indices:
                    try:
                        # DepMap values are Log2(TPM+1)
                        val = float(parts[idx])
                    except (ValueError, IndexError):
                        val = 0.0
                    vals.append(val)

                # Create DataFrame
                # Gene: The cleaned symbol
                # expression_norm: The value from the file (Log2(TPM+1))
                # expression_TPM: Calculated Linear TPM (2^x - 1)

                df = pd.DataFrame({
                    'Gene': clean_genes,
                    'expression_norm': vals
                })

                # Calculate linear TPM
                df['expression_TPM'] = (2 ** df['expression_norm']) - 1

                # Sort by expression (optional, but nice for inspection)
                df = df.sort_values('expression_norm', ascending=False)

                # Save
                out_name = f"{curr_id}_expression.csv"
                df.to_csv(os.path.join(output_dir, out_name), index=False)

                found_count += 1
                if found_count == len(target_ids):
                    break

    click.echo(f"✓ Processed {found_count} cell lines. Data in {output_dir}")


@main.command()
@click.option('--top-n', default=300, help='Genes for specific cell lines (Default: 300).')
@click.option('--top-n-generic', default=500, help='Genes per cell for Generic pool (Default: 500).')
@click.option('--genome', default='GRCh38', help='Genome version.')
@click.option('--force', is_flag=True, help='Overwrite existing cai_weights.json if it exists.')
def build_cai_weights(top_n, top_n_generic, genome, force):
    """
    Generates CAI weight profiles for the cohort and a Generic fallback.
    """
    data_dir = get_data_dir()
    out_path = os.path.join(data_dir, "cai_weights.json")

    if os.path.exists(out_path) and not force:
        click.echo(click.style(f"✓ CAI weights already exist at {out_path}. Use --force to recalculate.", fg="green"))
        return

    manifest_path = os.path.join(data_dir, "cell_cohort.json")
    expression_dir = os.path.join(data_dir, "processed_expression")

    if not os.path.exists(manifest_path):
        click.echo(click.style("❌ Error: cell_cohort.json not found.", fg="red"))
        return

    with open(manifest_path, 'r') as f:
        cohort = json.load(f)

    paths = get_paths(genome)
    # 1. Initialize Mapper first to get valid gene set
    mapper = GeneCoordinateMapper(paths['db'])
    valid_db_genes = set(mapper.gene_name_map.keys())

    # 2. PHASE 1: Use the fixed orchestrator function
    # Note: Pass valid_db_genes here!

    cell_line_top_genes, fallback_genes, global_fetch_set = load_cell_line_gene_maps(
        cell_map=cohort,
        data_dir=Path(expression_dir),
        valid_db_genes=valid_db_genes,
        n_specific=top_n_generic,
        n_fallback_scan=top_n_generic
    )

    # --- PHASE 2: Build Sequence Registry ---
    all_target_genes = set().union(*cell_line_top_genes.values(), global_fetch_set)
    click.echo(f"Fetching sequences for {len(all_target_genes)} valid genes...")

    gene_to_data = get_locus_to_data_dict(include_introns=False, gene_subset=list(all_target_genes))

    # Reuse the mapper we created above to save memory/time
    ref_registry = build_gene_sequence_registry(
        genes=list(all_target_genes),
        gene_to_data=gene_to_data,
        mapper=mapper
    )

    # --- PHASE 3: Calculate Weights ---
    cai_weights_map = {}

    # A. Specific Weights
    click.echo("\nCalculating cell-specific weights (Top 300)...")
    for cell_name, genes in cell_line_top_genes.items():
        cds_list = [
            ref_registry[g]['cds_sequence'] for g in genes
            if g in ref_registry and ref_registry[g].get('cds_sequence')
        ]
        cds_list = cds_list[:top_n]


        if len(cds_list) == top_n:
            _, weights_flat = calc_CAI_weight(cds_list)
            cai_weights_map[cell_name] = weights_flat
            click.echo(f"  ✓ {cell_name} weights ready ({len(cds_list)} genes).")
        else:
            click.echo(f"  ⚠ {cell_name} has insufficient sequences.")

    # B. Generic Weights (Consensus Intersection)
    click.echo(f"\nCalculating Generic Weights (Intersection of Top {top_n_generic})...")

    # Extract sequences for the consensus genes
    fallback_cds_list = [
        ref_registry[g]['cds_sequence'] for g in fallback_genes
        if g in ref_registry and ref_registry[g].get('cds_sequence')
    ]

    fallback_cds_list = fallback_cds_list[:top_n]

    if fallback_cds_list:
        _, generic_weights = calc_CAI_weight(fallback_cds_list)
        cai_weights_map['Generic'] = generic_weights
        click.echo(f"  ✓ Generic weights ready ({len(fallback_cds_list)} genes).")
    else:
        click.echo("  ⚠ Warning: No sequences found for fallback genes.")

    # --- PHASE 4: Save ---
    with open(out_path, 'w') as f:
        json.dump(cai_weights_map, f, indent=4)

    click.echo(click.style(f"\nSUCCESS: CAI weights saved to {out_path}", fg="green"))


@main.command()
@click.option('--force', is_flag=True, help="Force redownload if file exists.")
def setup_mrna_halflife(force):
    """
    Downloads the 'species_stability_no_threshold.csv.gz' dataset from the TTDB source.
    This file contains mRNA half-life data required for the stability features.
    """
    # The specific File ID
    FILE_ID = "1GekvDui-B2tSAQ6wgO3tIXpKd54EGRbn"

    # gdown works best with the 'uc' (User Content) export URL format
    url = f'https://drive.google.com/uc?id={FILE_ID}'

    # Assuming get_data_dir() is defined in your utils or imported
    data_dir = get_data_dir()

    # Ensure the directory exists
    os.makedirs(data_dir, exist_ok=True)

    destination = os.path.join(data_dir, 'mrna_half_life.csv.gz')

    click.echo(f"Initializing Stability Data setup...")
    click.echo(f"Target path: {destination}")

    if os.path.exists(destination) and not force:
        click.echo(click.style(f"✓ File already exists at {destination}. Use --force to overwrite.", fg="green"))
        return

    try:
        click.echo("Contacting Google Drive via gdown...")

        # gdown.download automatically handles the "virus scan warning" confirmation
        # quiet=False allows you to see the progress bar
        output = gdown.download(url, destination, quiet=False)

        if not output:
            raise Exception("Download failed (no output file).")

        click.echo(click.style(f"✓ Download complete: {destination}", fg="green"))

        # Verify it is a valid gzip (Crucial check for partial downloads)
        try:
            with gzip.open(destination, 'rb') as f:
                f.read(1)
            click.echo(click.style(f"✓ Integrity check passed (valid gzip).", fg="green"))
        except Exception:
            click.echo(
                click.style(f"⚠ Warning: The downloaded file is not a valid gzip. It might be an HTML error page.",
                            fg="red"))
            click.echo(click.style(f"Try deleting {destination} and running with --force.", fg="yellow"))
            sys.exit(1)

    except Exception as e:
        click.echo(click.style(f"Error downloading file: {e}", fg="red"))
        sys.exit(1)


GENCODE_HUMAN_RELEASE = "38"


def get_genome_metadata(genome):
    """
    Returns download URLs for supported genomes.
    """
    CONFIG = {
        # --- MAMMALS (GENCODE) ---
        'GRCh38': {
            'base_url': f'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{GENCODE_HUMAN_RELEASE}',
            'fasta_name' : 'GRCh38.p13.genome.fa.gz',
            'gtf_name': f'gencode.v{GENCODE_HUMAN_RELEASE}.chr_patch_hapl_scaff.annotation.gtf.gz'
        },
        'GRCm39': {
            'base_url': 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33',
            'fasta_name': 'GRCm39.genome.fa.gz',
            'gtf_name': 'gencode.vM33.chr_patch_hapl_scaff.basic.annotation.gtf.gz'
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