import gzip
import hashlib
import io
import itertools
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import time
from importlib.resources import files
from pathlib import Path

import click
import gdown
import gffutils
import pandas as pd
import requests
from gffutils.iterators import DataIterator
from pyfaidx import Fasta

from tauso.data.data import get_data_dir, get_paths
from tauso.features.codon_usage.cai import calc_CAI_weight
from tauso.features.codon_usage.find_cai_reference import load_cell_line_gene_maps
from tauso.genome.read_human_genome import get_locus_to_data_dict
from tauso.genome.TranscriptMapper import (
    GeneCoordinateMapper,
    build_gene_sequence_registry,
)
from tauso.off_target.search import find_all_gene_off_targets, get_bowtie_index_base

logger = logging.getLogger(__name__)


@click.group()
def main():
    """Tauso: ASO Design Toolkit"""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )


def download_and_gunzip(url, dest_path, remove_gz=False):
    if os.path.exists(dest_path):
        click.echo(f"  File already exists: {os.path.basename(dest_path)}")
        return

    try:
        click.echo(f"  Downloading {os.path.basename(url)}...")
        temp_gz = dest_path + ".gz"

        with requests.get(url, stream=True, timeout=30) as r:
            r.raise_for_status()
            total_size = int(r.headers.get("content-length", 0))

            with open(temp_gz, "wb") as f:
                with click.progressbar(length=total_size, label="    Downloading") as bar:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
                        bar.update(len(chunk))

        click.echo(f"  Unzipping to {os.path.basename(dest_path)}...")
        with gzip.open(temp_gz, "rb") as f_in:
            with open(dest_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        if remove_gz:
            os.remove(temp_gz)
    except Exception as e:
        if os.path.exists(dest_path):
            os.remove(dest_path)
        if os.path.exists(temp_gz):
            os.remove(temp_gz)
        raise e


def count_lines(filepath):
    with open(filepath, "rb") as f:
        return sum(1 for line in f if not line.startswith(b"#"))


def batch_iterator(iterator, batch_size=1000):
    while True:
        batch = list(itertools.islice(iterator, batch_size))
        if not batch:
            break
        yield batch


@main.command()
@click.option("--force", is_flag=True, help="Force redownload.")
def setup_depmap(force):
    """
    Downloads DepMap Public 25Q3 data directly from the DepMap API.
    Auto-fetches fresh signed URLs to ensure downloads work.
    OmicsExpression is converted to Parquet after download; original CSV is deleted.
    """
    data_dir = get_data_dir()
    os.makedirs(data_dir, exist_ok=True)

    RELEASE = "DepMap Public 25Q3"
    TARGET_FILES = {
        "Model.csv": "Model.csv",
        "OmicsProfiles.csv": "OmicsProfiles.csv",
        "OmicsExpressionTPMLogp1HumanAllGenesStranded.csv": "OmicsExpressionTPMLogp1HumanAllGenesStranded.csv",
    }
    # These CSVs are converted to Parquet and the original CSV deleted afterwards.
    CONVERT_TO_PARQUET = {"OmicsExpressionTPMLogp1HumanAllGenesStranded.csv"}

    click.echo(f"Initializing DepMap setup for: {RELEASE}")

    # Check what needs downloading before hitting the API
    to_download = {
        remote: local
        for remote, local in TARGET_FILES.items()
        if force or not os.path.exists(os.path.join(data_dir, local))
    }

    if not to_download:
        click.echo("✓ All DepMap files already present.")
        click.echo("\nDepMap setup complete.")
        return

    # Only fetch the index if we actually need to download something
    click.echo("Fetching fresh download URLs from DepMap API...")

    headers = {"User-Agent": "Mozilla/5.0"}
    index_url = "https://depmap.org/portal/api/download/files"

    try:
        r = requests.get(index_url, headers=headers)
        r.raise_for_status()
        index_df = pd.read_csv(io.StringIO(r.text))
    except Exception as e:
        click.echo(click.style(f"❌ Error fetching DepMap index: {e}", fg="red"))
        sys.exit(1)

    release_df = index_df[index_df["release"] == RELEASE]
    if release_df.empty:
        click.echo(click.style(f"❌ Release '{RELEASE}' not found in API.", fg="red"))
        click.echo(f"Available releases: {index_df['release'].unique()[:5]}...")
        sys.exit(1)

    for remote_name, local_name in TARGET_FILES.items():
        dest = os.path.join(data_dir, local_name)

        if local_name in CONVERT_TO_PARQUET:
            parquet_dest = dest.replace(".csv", ".parquet")
            # Migration: CSV exists but parquet does not → convert and delete
            if os.path.exists(dest) and not os.path.exists(parquet_dest):
                click.echo(f"  Converting {local_name} to Parquet...")
                pd.read_csv(dest).to_parquet(parquet_dest, index=False)
                os.remove(dest)
                click.echo(click.style(f"✓ Converted {local_name} → parquet, original deleted.", fg="green"))
            if os.path.exists(parquet_dest) and not force:
                click.echo(f"✓ {local_name} (parquet) exists.")
                continue
        else:
            if os.path.exists(dest) and not force:
                click.echo(f"✓ {local_name} exists.")
                continue

        file_row = release_df[release_df["filename"] == remote_name]
        if file_row.empty:
            click.echo(click.style(f"⚠ Warning: '{remote_name}' not found in {RELEASE}.", fg="yellow"))
            continue

        download_url = file_row.iloc[0]["url"]

        click.echo(f"Downloading {local_name}...")
        try:
            subprocess.run(
                ["wget", "-q", "--show-progress", "-U", "Mozilla/5.0", "-O", dest, download_url],
                check=True,
            )
        except subprocess.CalledProcessError:
            click.echo(click.style(f"❌ Failed to download {local_name}", fg="red"))
            if os.path.exists(dest):
                os.remove(dest)
            continue

        click.echo(click.style(f"✓ Downloaded {local_name}.", fg="green"))

        if local_name in CONVERT_TO_PARQUET:
            click.echo(f"  Converting {local_name} to Parquet...")
            pd.read_csv(dest).to_parquet(parquet_dest, index=False)
            os.remove(dest)
            click.echo(click.style(f"✓ Converted to Parquet, original CSV deleted.", fg="green"))

    click.echo("\nDepMap setup complete.")


import re


@main.command()
@click.argument("cell_names", nargs=-1)
@click.option("--reset", is_flag=True, help="Clear existing list before adding.")
def add_cell(cell_names, reset):
    """
    Search for cell lines by name and add them to the project cohort.
    Example: tauso add-cell "HepG2" "HeLa" "U251"
    """
    data_dir = get_data_dir()
    profiles_path = os.path.join(data_dir, "OmicsProfiles.csv")
    manifest_path = os.path.join(data_dir, "cell_cohort.json")

    if not os.path.exists(profiles_path):
        click.echo(
            click.style(
                "Error: OmicsProfiles.csv not found. Run 'setup-depmap' first.",
                fg="red",
            )
        )
        sys.exit(1)

    # 1. Load Existing Manifest
    cohort = {}
    if os.path.exists(manifest_path) and not reset:
        with open(manifest_path, "r") as f:
            cohort = json.load(f)

    # 2. Load Metadata (Optimized)
    click.echo("Loading metadata...")
    df = pd.read_csv(
        profiles_path,
        usecols=["ModelID", "StrippedCellLineName", "IsDefaultEntryForModel"],
    )
    # Normalize for fuzzy search
    df["clean"] = df["StrippedCellLineName"].str.replace(r"[^a-zA-Z0-9]", "", regex=True).str.upper()

    # 3. Search
    for query in cell_names:
        # FIXED: Use the exact same regex as the dataframe to avoid mismatch on periods/underscores
        clean_q = re.sub(r"[^a-zA-Z0-9]", "", query).upper()

        # Exact match
        match = df[df["clean"] == clean_q]
        # Fuzzy match
        if match.empty:
            match = df[df["clean"].str.contains(clean_q, na=False)]

        if not match.empty:
            # Prefer default entry
            best = match[match["IsDefaultEntryForModel"] == True]
            if best.empty:
                best = match.iloc[[0]]

            found_name = best.iloc[0]["StrippedCellLineName"]
            ach_id = best.iloc[0]["ModelID"]

            cohort[found_name] = ach_id
            click.echo(click.style(f"✓ Found: {query} -> {found_name} ({ach_id})", fg="green"))
        else:
            # ADDED: Reasons why it might not be found
            reasons = " (Reasons: Known by different alias, primary/non-cancer cell line, or lacks Omics profiling)"
            click.echo(click.style(f"⚠ Not Found: {query}{reasons}", fg="yellow"))

    # 4. Save
    with open(manifest_path, "w") as f:
        json.dump(cohort, f, indent=4)

    click.echo(f"Cohort saved to {manifest_path} ({len(cohort)} cell lines).")


@main.command()
@click.option("--genome", default="GRCh38", help="Genome version (default: GRCh38).")
def build_omics(genome):
    """
    Generates gene-level expression files for all cell lines in the cohort.
    Uses DepMap 'AllGenes' format (Gene Level Log2(TPM+1)).
    """
    if genome != "GRCh38":
        # While the logic is genome-agnostic now, we keep this guard if your downstream tools expect GRCh38
        click.echo(
            click.style(
                f"⚠ Warning: This command is optimized for DepMap (Human GRCh38) data.",
                fg="yellow",
            )
        )

    data_dir = get_data_dir()
    manifest_path = os.path.join(data_dir, "cell_cohort.json")
    exp_path = os.path.join(data_dir, "OmicsExpressionTPMLogp1HumanAllGenesStranded.parquet")

    if not os.path.exists(manifest_path):
        click.echo(click.style("❌ No cohort found. Use 'tauso add-cell' first.", fg="red"))
        return

    if not os.path.exists(exp_path):
        click.echo(click.style(f"❌ Expression parquet not found: {exp_path}", fg="red"))
        click.echo("Run 'tauso setup-depmap' to download and convert it.")
        return

    with open(manifest_path, "r") as f:
        cohort = json.load(f)

    target_ids = set(cohort.values())
    click.echo(f"Processing {len(target_ids)} cell lines from cohort...")

    click.echo(f"Loading {os.path.basename(exp_path)}...")
    exp_df = pd.read_parquet(exp_path)

    model_col = "ModelID" if "ModelID" in exp_df.columns else exp_df.columns[0]
    gene_cols = [c for c in exp_df.columns if c != model_col]

    gene_regex = re.compile(r"^(.+?) \(\d+\)$")
    clean_gene_map = {c: (m.group(1) if (m := gene_regex.match(c)) else c) for c in gene_cols}

    output_dir = os.path.join(data_dir, "processed_expression")
    os.makedirs(output_dir, exist_ok=True)

    found_count = 0
    for curr_id in target_ids:
        cell_rows = exp_df[exp_df[model_col] == curr_id]
        if cell_rows.empty:
            continue

        click.echo(f"  Extracting {curr_id}...")
        row = cell_rows.iloc[0]
        vals = pd.to_numeric(row[gene_cols], errors="coerce").fillna(0.0).values
        clean_genes = [clean_gene_map[c] for c in gene_cols]

        out_df = pd.DataFrame({"Gene": clean_genes, "expression_norm": vals})
        out_df["expression_TPM"] = (2 ** out_df["expression_norm"]) - 1
        out_df = out_df.sort_values("expression_norm", ascending=False)
        out_df.to_csv(os.path.join(output_dir, f"{curr_id}_expression.csv"), index=False)
        found_count += 1

    click.echo(f"✓ Processed {found_count} cell lines. Data in {output_dir}")


@main.command()
@click.option("--top-n", default=300, help="Genes for specific cell lines (Default: 300).")
@click.option(
    "--top-n-generic",
    default=500,
    help="Genes per cell for Generic pool (Default: 500).",
)
@click.option("--genome", default="GRCh38", help="Genome version.")
@click.option("--force", is_flag=True, help="Overwrite existing cai_weights.json if it exists.")
def build_cai_weights(top_n, top_n_generic, genome, force):
    """
    Generates CAI weight profiles for the cohort and a Generic fallback.
    """
    data_dir = get_data_dir()
    out_path = os.path.join(data_dir, "cai_weights.json")

    if os.path.exists(out_path) and not force:
        click.echo(
            click.style(
                f"✓ CAI weights already exist at {out_path}. Use --force to recalculate.",
                fg="green",
            )
        )
        return

    manifest_path = os.path.join(data_dir, "cell_cohort.json")
    expression_dir = os.path.join(data_dir, "processed_expression")

    if not os.path.exists(manifest_path):
        click.echo(click.style("❌ Error: cell_cohort.json not found.", fg="red"))
        return

    with open(manifest_path, "r") as f:
        cohort = json.load(f)

    paths = get_paths(genome)
    # 1. Initialize Mapper first to get valid gene set
    mapper = GeneCoordinateMapper(paths["gtf_db"])
    valid_db_genes = set(mapper.gene_name_map.keys())

    # 2. PHASE 1: Use the fixed orchestrator function
    # Note: Pass valid_db_genes here!

    cell_line_top_genes, fallback_genes, global_fetch_set = load_cell_line_gene_maps(
        cell_map=cohort,
        data_dir=Path(expression_dir),
        valid_db_genes=valid_db_genes,
        n_specific=top_n_generic,
        n_fallback_scan=top_n_generic,
        filter_mode="protein_coding",
        genome_db=mapper.db,
    )

    # --- PHASE 2: Build Sequence Registry ---
    all_target_genes = set().union(*cell_line_top_genes.values(), global_fetch_set)
    click.echo(f"Fetching sequences for {len(all_target_genes)} valid genes...")

    gene_to_data = get_locus_to_data_dict(include_introns=False, gene_subset=list(all_target_genes))

    # Reuse the mapper we created above to save memory/time
    ref_registry = build_gene_sequence_registry(genes=list(all_target_genes), gene_to_data=gene_to_data, mapper=mapper)

    # --- PHASE 3: Calculate Weights ---
    cai_weights_map = {}

    # A. Specific Weights
    click.echo("\nCalculating cell-specific weights (Top 300)...")
    for cell_name, genes in cell_line_top_genes.items():
        cds_list = [
            ref_registry[g]["cds_sequence"] for g in genes if g in ref_registry and ref_registry[g].get("cds_sequence")
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
        ref_registry[g]["cds_sequence"]
        for g in fallback_genes
        if g in ref_registry and ref_registry[g].get("cds_sequence")
    ]

    fallback_cds_list = fallback_cds_list[:top_n]

    if fallback_cds_list:
        _, generic_weights = calc_CAI_weight(fallback_cds_list)
        cai_weights_map["Generic"] = generic_weights
        click.echo(f"  ✓ Generic weights ready ({len(fallback_cds_list)} genes).")
    else:
        click.echo("  ⚠ Warning: No sequences found for fallback genes.")

    # --- PHASE 4: Save ---
    with open(out_path, "w") as f:
        json.dump(cai_weights_map, f, indent=4)

    click.echo(click.style(f"\nSUCCESS: CAI weights saved to {out_path}", fg="green"))


def calc_file_hash(fpath):
    sha256_hash = hashlib.sha256()
    with open(fpath, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()


def should_skip_download(file_path, force, expected_hash):
    """Checks if file exists and matches hash. Returns True if we can skip."""
    if not os.path.exists(file_path):
        return False

    if force:
        click.echo(click.style("Force flag detected. Redownloading...", fg="yellow"))
        return False

    click.echo("File exists. Verifying hash...")
    try:
        current_hash = calc_file_hash(file_path)
        if current_hash == expected_hash:
            click.echo(click.style(f"✓ Hash matched. Skipping download.", fg="green"))
            return True
        else:
            click.echo(
                click.style(
                    f"⚠ Hash mismatch on existing file. Redownload with --force",
                    fg="red",
                )
            )
            click.echo(f"Expected: {expected_hash}")
            click.echo(f"Got:      {current_hash}")
            sys.exit(1)
    except Exception as e:
        click.echo(
            click.style(
                f"⚠ Error reading existing file ({e}). Redownload with --force",
                fg="red",
            )
        )
        sys.exit(1)


@main.command()
@click.option("--force", is_flag=True, help="Force redownload if file exists.")
def setup_mrna_halflife(force):
    """
    Downloads the 'species_stability_no_threshold.csv.gz' dataset from the TTDB source.
    """
    FILE_ID = "1GekvDui-B2tSAQ6wgO3tIXpKd54EGRbn"
    EXPECTED_SHA256 = "ec0c1f90eed96516de510c484a7a7dd4bbd10d253306710bb33e714f88c5c135"

    url = f"https://drive.google.com/uc?id={FILE_ID}"
    data_dir = get_data_dir()
    os.makedirs(data_dir, exist_ok=True)
    destination = os.path.join(data_dir, "mrna_half_life.csv.gz")

    click.echo(f"Initializing Stability Data setup...")
    click.echo(f"Target path: {destination}")

    # 1. Check if we can skip
    if should_skip_download(destination, force, EXPECTED_SHA256):
        return

    # 2. Perform Download
    try:
        click.echo("Contacting Google Drive via gdown...")
        output = gdown.download(url, destination, quiet=False)

        if not output:
            raise Exception("Download failed (no output file).")

        click.echo(click.style(f"✓ Download complete: {destination}", fg="green"))

        # 3. Verify Gzip Integrity
        try:
            with gzip.open(destination, "rb") as f:
                f.read(1)
            click.echo(click.style(f"✓ Integrity check passed (valid gzip).", fg="green"))
        except Exception:
            click.echo(
                click.style(
                    f"⚠ Warning: File is not a valid gzip (likely HTML error).",
                    fg="red",
                )
            )
            sys.exit(1)

        # 4. Verify Hash (Post-Download)
        click.echo("Verifying new file hash...")
        new_hash = calc_file_hash(destination)

        if new_hash != EXPECTED_SHA256:
            click.echo(click.style(f"⚠ Hash mismatch!", fg="red"))
            click.echo(f"Expected: {EXPECTED_SHA256}")
            click.echo(f"Got:      {new_hash}")
            click.echo(
                click.style(
                    "Please update EXPECTED_SHA256 with the 'Got' value above.",
                    fg="yellow",
                )
            )
            sys.exit(1)

        click.echo(click.style(f"✓ Hash check passed.", fg="green"))

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
        "GRCh38": {
            "base_url": f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{GENCODE_HUMAN_RELEASE}",
            "fasta_name": "GRCh38.primary_assembly.genome.fa.gz",
            "gtf_name": f"gencode.v{GENCODE_HUMAN_RELEASE}.annotation.gtf.gz",
            "gff_name": f"gencode.v{GENCODE_HUMAN_RELEASE}.annotation.gff3.gz",
        },
        "GRCm39": {
            "base_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33",
            "fasta_name": "GRCm39.primary_assembly.genome.fa.gz",
            "gtf_name": "gencode.vM33.annotation.gtf.gz",
        },
        # --- YEAST (Ensembl) ---
        "R64-1-1": {
            "base_url": "http://ftp.ensembl.org/pub/release-110",  # Stable release
            "fasta_name": "fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz",
            "gtf_name": "gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.110.gtf.gz",
        },
        # --- E. COLI (Ensembl Bacteria) ---
        # Note: E. coli K-12 MG1655
        "ASM584v2": {
            "base_url": "http://ftp.ensemblgenomes.org/pub/bacteria/release-57",
            "fasta_name": "fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.toplevel.fa.gz",
            "gtf_name": "gtf/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.57.gtf.gz",
        },
    }

    if genome not in CONFIG:
        return None

    meta = CONFIG[genome]
    base = meta["base_url"]

    url_dict = dict()

    if "fasta_name" in meta:
        url_dict["fasta"] = f"{base}/{meta['fasta_name']}"
    if "gtf_name" in meta:
        url_dict["gtf"] = f"{base}/{meta['gtf_name']}"
    if "gff_name" in meta:
        url_dict["gff"] = f"{base}/{meta['gff_name']}"

    return url_dict


def build_annotation_db(
    annotation_path: str,
    db_path: str,
    db_success_path: str,
    genome: str,
    fmt: str = "GTF",  # "GTF" or "GFF"
) -> None:
    """
    Build a gffutils database from a GTF or GFF annotation file.

    Args:
        annotation_path:  Path to the source .gtf / .gff file.
        db_path:          Path where the SQLite database will be written.
        db_success_path:  Path to the sentinel file written on success.
        genome:           Genome label used in user-facing messages.
        fmt:              Format label shown in progress messages ("GTF" or "GFF").

    Raises:
        SystemExit: on any database-build failure (after cleanup).
    """
    # ── Guard: already built ────────────────────────────────────────────────
    if os.path.exists(db_path):
        if os.path.exists(db_success_path):
            click.echo(f"✓ {fmt} Database already exists at {db_path}")
            return
        click.echo(
            click.style(
                f"⚠ Found incomplete/corrupt {fmt} database. Rebuilding...",
                fg="yellow",
            )
        )
        os.remove(db_path)

    # ── Build ────────────────────────────────────────────────────────────────
    try:
        click.echo(f"Building database at {db_path}...")

        total_lines = count_lines(annotation_path)
        data_it = DataIterator(annotation_path)

        click.echo(f"  - Parsing {total_lines:,} annotation lines...")
        with click.progressbar(length=total_lines, label=f"    Parsing {fmt}") as bar:

            def progress_wrapper(iterator=data_it, progress_bar=bar):
                for feature in iterator:
                    progress_bar.update(1)
                    yield feature

            merge_strat = "merge" if fmt == "GTF" else "create_unique"
            keep_ord = fmt == "GTF"
            sort_attrs = fmt == "GTF"

            db = gffutils.create_db(
                progress_wrapper(),
                dbfn=db_path,
                force=True,
                keep_order=keep_ord,
                merge_strategy=merge_strat,
                sort_attribute_values=sort_attrs,
                disable_infer_genes=True,
                disable_infer_transcripts=True,
            )

        print()  # newline after progress bar

        # ── Write success sentinel ───────────────────────────────────────────
        with open(db_success_path, "w") as f:  # ← fixed: was writing to db_path
            f.write("Setup completed successfully.")

        click.echo(click.style(f"✓ Setup complete for {genome}.", fg="green"))

    except Exception as e:
        click.echo(click.style(f"\nError building {fmt} database: {e}", fg="red"))
        for path in (db_path, db_success_path):
            if os.path.exists(path):
                os.remove(path)
        sys.exit(1)


@main.command()
@click.option("--genome", default="GRCh38", help="Genome name (GRCh38 or GRCm39).")
@click.option("--force", is_flag=True, help="Force re-download and rebuild.")
@click.option("--remove-gz", is_flag=True, help="Remove the zipped files after download.")
def setup_genome(genome, force, remove_gz):
    """
    Sets up the genome environment (Download -> Index -> Database).
    Supports GRCh38 (Human) and GRCm39 (Mouse).
    Always downloads the 'basic' gene annotation subset.

    :Note: You may set the TAUSO_DATA_DIR environment variable to switch the folder.
    """
    paths = get_paths(genome)
    fasta_path = paths["fasta"]
    gtf_path = paths["gtf"]
    gff_path = paths["gff"]
    gtf_db_path = paths["gtf_db"]
    gff_db_path = paths["gff_db"]

    gtf_db_success = gtf_db_path + ".success"
    gff_db_success = gff_db_path + ".success"

    if force:
        click.echo("Force flag detected. Cleaning up old files...")
        for f in [
            fasta_path,
            gtf_path,
            gff_path,
            gtf_db_path,
            gff_db_path,
            fasta_path + ".fai",
            gtf_db_success,
            gff_db_success,
        ]:
            if os.path.exists(f):
                os.remove(f)

    # --- PHASE 1: Download & Index ---
    try:
        # Check if we have automatic download support for this genome
        url_dict = get_genome_metadata(genome)

        if url_dict:
            for file_type, download_url in url_dict.items():
                if file_type == "gff" and genome != "GRCh38":  # TODO: support GFF for all
                    continue
                db_key = f"{file_type}_db"
                if db_key in paths:
                    db_success = paths[db_key] + ".success"
                    if os.path.exists(paths[db_key]) and os.path.exists(db_success):
                        click.echo(f"✓ {file_type} database already built, skipping download.")
                        continue
                if not os.path.exists(paths[file_type]):
                    # Check for a cached compressed copy before hitting the network
                    gz_key = f"{file_type}_gz"
                    gz_path = paths.get(gz_key, paths[file_type] + ".gz")
                    if os.path.exists(gz_path):
                        click.echo(f"  Found cached {os.path.basename(gz_path)}, decompressing...")
                        with gzip.open(gz_path, "rb") as f_in, open(paths[file_type], "wb") as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    else:
                        click.echo(f"Downloading {file_type} from GENCODE...")
                        download_and_gunzip(download_url, paths[file_type], remove_gz=remove_gz)
        else:
            # Fallback for unsupported genomes
            if not os.path.exists(fasta_path) or not os.path.exists(gtf_path):
                click.echo(
                    click.style(
                        f"Error: Automatic download not supported for '{genome}'.",
                        fg="red",
                    )
                )
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

    # --- PHASE 3: Build GTF & GFF Databases ---
    # (Database building logic remains identical)
    build_annotation_db(
        annotation_path=gtf_path, db_path=gtf_db_path, db_success_path=gtf_db_success, genome=genome, fmt="GTF"
    )
    if genome == "GRCh38":
        build_annotation_db(
            annotation_path=gff_path, db_path=gff_db_path, db_success_path=gff_db_success, genome=genome, fmt="GFF"
        )


# --- NEW COMMAND: OFF-TARGET SEARCH ---
@main.command()
@click.argument("sequence")
@click.option("--genome", default="GRCh38", help="Genome name (default: GRCh38).")
@click.option("--mismatches", "-m", default=3, help="Max mismatches allowed.")
@click.option("--output", "-o", default=None, help="Output CSV file.")
def run_off_target(sequence, genome, mismatches, output):
    """
    Search for off-target binding sites for a given ASO sequence.
    """
    try:
        paths = get_paths(genome)
        # Construct the expected path to the "SUCCESS" sentinel file
        # This mirrors the logic in get_bowtie_index_base
        fasta_dir = os.path.dirname(paths["fasta"])
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
    hits_df = hits_df.sort_values(by=["mismatches", "chrom", "start"])

    # Display Top 20
    click.echo("\nTop Hits:")
    click.echo(f"{'Chrom':<10} | {'Start':<12} | {'Str':<3} | {'MM':<2} | {'Gene':<15} | {'Region'}")
    click.echo("-" * 70)

    count = 0
    for _, hit in hits_df.iterrows():
        if count >= 20:
            click.echo(f"... and {len(hits_df) - 20} more.")
            break

        gene = str(hit["gene_name"]) if hit["gene_name"] else "None"
        click.echo(
            f"{hit['chrom']:<10} | {hit['start']:<12} | {hit['strand']:<3} | {hit['mismatches']:<2} | {gene:<15} | {hit['region_type']}"
        )
        count += 1

    if output:
        hits_df.to_csv(output, index=False)
        click.echo(f"\nFull results saved to: {output}")


@main.command()
@click.option("--genome", default="GRCh38", help="Genome name (default: GRCh38).")
@click.option("--force", is_flag=True, help="Force rebuild of the index.")
@click.option("--threads", "-t", default=1, help="Number of CPU threads to use (default: 1).")
@click.option("--mem-per-thread", default=800, help="Max memory limit per thread in MB (default: 800).")
def setup_bowtie(genome, force, threads, mem_per_thread):
    """
    Generates (or validates) the Bowtie index for the specified genome.
    Requires 'setup-genome' to be run first to ensure the FASTA file exists.
    """
    click.echo(f"Initializing Bowtie setup for {genome} (Threads: {threads}, Max Mem/Thread: {mem_per_thread}MB)...")

    try:
        index_path = get_bowtie_index_base(
            genome=genome, force_rebuild=force, threads=threads, mem_per_thread_mb=mem_per_thread
        )

        click.echo(click.style(f"✓ Bowtie index ready at: {index_path}", fg="green"))

    except FileNotFoundError:
        click.echo(click.style(f"Error: FASTA file for {genome} not found.", fg="red"))
        click.echo(f"Please run 'tauso setup-genome --genome {genome}' first.")
        sys.exit(1)
    except Exception as e:
        click.echo(click.style(f"Error building Bowtie index: {e}", fg="red"))
        sys.exit(1)


@main.command(context_settings=dict(ignore_unknown_options=True, allow_extra_args=True))
@click.pass_context
@click.option(
    "--force-clone",
    "-f",
    is_flag=True,
    help="Force re-cloning of the raccess repository (passed to install_raccess.sh).",
)
def setup_raccess(ctx, force_clone):
    """
    Run the raccess installation per their license.
    Installs into the configured TAUSO_DATA_DIR.
    """
    script_path = files("tauso._raccess") / "install_raccess.sh"
    data_dir = get_data_dir()
    raccess_dir = os.path.join(data_dir, "raccess")

    forwarded_args = list(ctx.args)
    if force_clone:
        forwarded_args.append("--force-clone")

    # Pass raccess_dir as the first argument to the script
    cmd = ["bash", str(script_path), raccess_dir] + forwarded_args

    click.echo(f"+ {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logger.error("install_raccess.sh failed, consider installing zlib1g-dev")
        sys.exit(1)


if __name__ == "__main__":
    main()
