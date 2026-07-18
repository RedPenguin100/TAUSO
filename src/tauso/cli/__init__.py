import gzip
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
import gffutils
import pandas as pd
from gffutils.iterators import DataIterator
from pyfaidx import Fasta

from tauso.cli_utils import (
    count_lines,
    download_and_gunzip,
    download_with_progress,
    echo_err,
    echo_ok,
    echo_warn,
    sha1_file,
    sha256_file,
    verify_hash_or_exit,
)
from tauso.data.data import get_data_dir, get_paths, load_gtf_db
from tauso.features.codon_usage.cai import CAI_DEFAULT_PSEUDOCOUNT, CAI_WEIGHTS_FILENAME, build_scorer_from_reference
from tauso.features.codon_usage.find_cai_reference import load_cell_line_gene_maps
from tauso.features.codon_usage.tai import TGCNSource
from tauso.features.context.mrna_halflife import HALFLIFE_SOURCE_COLUMNS
from tauso.genome.read_human_genome import build_locus_cache, get_locus_to_data_dict
from tauso.genome.TranscriptMapper import build_gene_sequence_registry
from tauso.off_target.search import find_all_gene_off_targets, get_bowtie_index_base

from ..util import rna_to_dna
from ._download import (
    DEPMAP_FILES_SHA1,
    RRNA_SHA1,
    ZENODO_RRNA_RECORD,
    _ensure_depmap_file,
    _ensure_zenodo_content_file,
    _ensure_zenodo_table,
    _zenodo_file_url,
)

logger = logging.getLogger(__name__)


@click.group()
def main():
    """Tauso: ASO Design Toolkit"""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )


@main.command()
@click.option("--force", is_flag=True, help="Force redownload.")
def setup_depmap(force):
    """
    Downloads DepMap Public 25Q3 data from a pinned Zenodo mirror
    (https://zenodo.org/records/20355477), verifying each file's SHA1.
    The OmicsExpression CSV is converted to Parquet for fast loading and the
    CSV is then removed — all consumers read the Parquet directly. A sidecar
    `<parquet>.sha256` is written on conversion; future runs skip the 1.1 GB
    CSV download iff the Parquet still matches that hash. A mismatch triggers
    a fresh download + re-conversion.
    """
    data_dir = get_data_dir()
    os.makedirs(data_dir, exist_ok=True)

    click.echo("Initializing DepMap setup (Zenodo mirror of DepMap Public 25Q3)...")

    omics_csv_name = "OmicsExpressionTPMLogp1HumanAllGenesStranded.csv"
    omics_csv = os.path.join(data_dir, omics_csv_name)
    omics_parquet = omics_csv.replace(".csv", ".parquet")
    omics_parquet_sha = omics_parquet + ".sha256"
    parquet_already_built = os.path.exists(omics_parquet) and not force

    # If we have both the Parquet and its sidecar hash, verify they agree.
    # Mismatch (or missing sidecar with --force) means the Parquet is no longer
    # the one we wrote → rebuild from a fresh CSV download.
    if parquet_already_built and os.path.exists(omics_parquet_sha):
        recorded = Path(omics_parquet_sha).read_text().strip()
        if sha256_file(omics_parquet) != recorded:
            echo_warn(f"{os.path.basename(omics_parquet)} hash mismatch — will re-download CSV and re-convert.")
            os.remove(omics_parquet)
            os.remove(omics_parquet_sha)
            parquet_already_built = False

    for filename, expected_sha1 in DEPMAP_FILES_SHA1.items():
        if filename == omics_csv_name and parquet_already_built:
            echo_ok(f"{filename} → Parquet already present; skipping CSV download.")
            continue
        _ensure_depmap_file(filename, expected_sha1, data_dir, force)

    if not parquet_already_built:
        click.echo("  Converting OmicsExpression CSV to Parquet...")
        pd.read_csv(omics_csv).to_parquet(omics_parquet, index=False)
        echo_ok("Converted to Parquet.")

    # Record (or refresh) the sidecar hash so future runs can detect tampering / bit rot.
    if not os.path.exists(omics_parquet_sha):
        Path(omics_parquet_sha).write_text(sha256_file(omics_parquet))

    # The CSV is dead weight once the Parquet exists — every consumer (production
    # code in genome/transcriptome.py and features/expression/general_expression.py,
    # plus tests) takes the CSV path but immediately swaps the suffix to .parquet.
    if os.path.exists(omics_csv):
        os.remove(omics_csv)
        echo_ok(f"Removed {omics_csv_name} (Parquet supersedes it).")

    click.echo("\nDepMap setup complete.")


@main.command(name="setup-omics")
@click.option("--force", is_flag=True, help="Force redownload of all omics datasets.")
@click.pass_context
def setup_omics(ctx, force):
    """
    Set up all omics datasets: DepMap (cell-line metadata, profiles, expression),
    mRNA half-life data, human tGCN (tAI), ATtRACT RBP motifs, and the ribo-seq
    bigWig. Use 'build-cohort-expression' afterwards to derive per-cohort
    expression files.
    """
    click.echo(click.style("=== setup-omics: DepMap ===", bold=True))
    ctx.invoke(setup_depmap, force=force)
    click.echo()
    click.echo(click.style("=== setup-omics: mRNA half-life ===", bold=True))
    ctx.invoke(setup_mrna_halflife, force=force)
    click.echo()
    click.echo(click.style("=== setup-omics: human tGCN ===", bold=True))
    ctx.invoke(setup_tgcn, force=force)
    click.echo()
    click.echo(click.style("=== setup-omics: ATtRACT RBP ===", bold=True))
    ctx.invoke(setup_attract, force=force)
    click.echo()
    click.echo(click.style("=== setup-omics: ribo-seq ===", bold=True))
    ctx.invoke(setup_riboseq, force=force)
    click.echo()
    echo_ok("Omics setup complete.")


@main.command(name="setup-all")
@click.option("--genome", default="GRCh38", help="Genome to set up (default: GRCh38).")
@click.option("--force", is_flag=True, help="Force redownload of every dataset.")
@click.option("--threads", "-t", default=1, help="Threads for bowtie index build (default: 1).")
@click.option("--mem-per-thread", default=800, help="Max MB/thread for bowtie (default: 800).")
@click.pass_context
def setup_all(ctx, genome, force, threads, mem_per_thread):
    """
    End-to-end setup: genome + bowtie + omics + raccess. Idempotent — already-present
    datasets are verified by hash and skipped unless --force is given.
    """
    click.echo(click.style("=== setup-all: genome ===", bold=True))
    ctx.invoke(setup_genome, genome=genome, force=force, remove_gz=False)
    click.echo()
    click.echo(click.style("=== setup-all: bowtie ===", bold=True))
    ctx.invoke(setup_bowtie, genome=genome, force=force, threads=threads, mem_per_thread=mem_per_thread)
    click.echo()
    click.echo(click.style("=== setup-all: omics ===", bold=True))
    ctx.invoke(setup_omics, force=force)
    click.echo()
    click.echo(click.style("=== setup-all: raccess ===", bold=True))
    ctx.invoke(setup_raccess)
    click.echo()
    click.echo(click.style("=== setup-all: rRNA ===", bold=True))
    ctx.invoke(setup_rrna, force=force)
    click.echo()
    echo_ok("setup-all complete.")


DEFAULT_COHORT_CELLS = (
    "HEPG2",
    "SNU449",
    "HELA",
    "A431",
    "SKMEL28",
    "SHSY5Y",
    "U251MG",
    "NCIH929",
    "KMS11",
    "NCIH460",
    "SKNAS",
    "SKNSH",
    "KARPAS299",
    "HEP3B217",
    "THP1",
    "LNCAPCLONEFGC",
    "T24",
    "A549",
    "VCAP",
    "HUH7",
    "JURKAT",
    "SKOV3",
    "K562",
    "A172",
    "PC3",
    "MCF7",
    "SW872",
    "G361",
    "HEK293",
    "MM1S",
)


@main.command(name="build-cell-context")
@click.argument("cell_names", nargs=-1)
@click.option("--reset", is_flag=True, help="Clear existing cohort before adding.")
@click.option("--genome", default="GRCh38", help="Genome version (default: GRCh38).")
@click.option("--cai/--no-cai", default=True, help="Build CAI weights (default: enabled).")
@click.option("--force", is_flag=True, help="Force CAI weights recomputation.")
@click.pass_context
def build_cell_context(ctx, cell_names, reset, genome, cai, force):
    """
    Build all per-cell-line reference data needed by downstream feature
    computation: register cells, extract per-cell expression files, and
    (by default) compute CAI weights.

    Bare invocation uses the default cohort. Pass cell names to use a
    custom cohort; pass --reset to replace an existing cohort; pass
    --no-cai to skip the CAI weight computation step.

    Example:
      tauso build-cell-context                            # default cohort, with CAI
      tauso build-cell-context --no-cai                   # skip CAI
      tauso build-cell-context HEPG2 SNU449 HELA          # custom cohort (append)
      tauso build-cell-context --reset HEPG2 SNU449       # custom cohort (replace)
    """
    data_dir = get_data_dir()
    cohort_exists = os.path.exists(os.path.join(data_dir, "cell_cohort.json"))

    if cell_names:
        cells_to_add = cell_names
    elif reset:
        echo_err("--reset requires at least one cell name.")
        sys.exit(1)
    elif not cohort_exists:
        click.echo(f"No existing cohort; using {len(DEFAULT_COHORT_CELLS)} default cells.")
        cells_to_add = DEFAULT_COHORT_CELLS
    else:
        cells_to_add = ()

    if cells_to_add:
        click.echo(click.style(f"=== build-cell-context: add-cell ({len(cells_to_add)} cells) ===", bold=True))
        ctx.invoke(add_cell, cell_names=cells_to_add, reset=reset)
        click.echo()

    click.echo(click.style("=== build-cell-context: cohort expression ===", bold=True))
    ctx.invoke(build_cohort_expression, genome=genome)
    click.echo()
    if cai:
        click.echo(click.style("=== build-cell-context: CAI weights ===", bold=True))
        ctx.invoke(build_cai_weights, top_n=300, top_n_generic=500, genome=genome, force=force)
        click.echo()

    click.echo(click.style("=== build-cell-context: human tGCN ===", bold=True))
    ctx.invoke(setup_tgcn, force=force)
    click.echo()
    echo_ok("build-cell-context complete.")


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
def build_cohort_expression(genome):
    """
    Generates gene-level expression files for all cell lines in the cohort.
    Uses DepMap 'AllGenes' format (Gene Level Log2(TPM+1)).
    """
    if genome != "GRCh38":
        # While the logic is genome-agnostic now, we keep this guard if your downstream tools expect GRCh38
        echo_warn("This command is optimized for DepMap (Human GRCh38) data.")

    data_dir = get_data_dir()
    manifest_path = os.path.join(data_dir, "cell_cohort.json")
    exp_path = os.path.join(data_dir, "OmicsExpressionTPMLogp1HumanAllGenesStranded.parquet")

    if not os.path.exists(manifest_path):
        echo_err("No cohort found. Use 'tauso add-cell' first.")
        return

    if not os.path.exists(exp_path):
        echo_err(f"Expression parquet not found: {exp_path}")
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
@click.option(
    "--pseudocount",
    default=CAI_DEFAULT_PSEUDOCOUNT,
    show_default=True,
    type=int,
    help="Pseudocount added to codon counts before normalisation (codonbias). "
    "0 matches the previous hand-rolled math; 1 (codonbias's own default) "
    "smooths and avoids zero weights for unobserved codons.",
)
@click.option("--force", is_flag=True, help="Overwrite existing cai_weights.json if it exists.")
def build_cai_weights(top_n, top_n_generic, genome, pseudocount, force):
    """
    Generates CAI weight profiles for the cohort and a Generic fallback,
    using codonbias.scores.CodonAdaptationIndex over the top-N highly
    expressed CDS sequences per cell line. Saves {cell_line: {codon: weight}}
    to cai_weights.json.
    """
    data_dir = get_data_dir()
    out_path = os.path.join(data_dir, CAI_WEIGHTS_FILENAME)

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

    # 1. Load the genome annotation DB and the set of valid gene names from it
    genome_db = load_gtf_db(genome)
    valid_db_genes = set()
    for gene in genome_db.features_of_type("gene"):
        names = gene.attributes.get("gene_name")
        if names:
            valid_db_genes.add(names[0])

    # 2. PHASE 1: Use the fixed orchestrator function
    cell_line_top_genes, fallback_genes, global_fetch_set = load_cell_line_gene_maps(
        cell_map=cohort,
        data_dir=Path(expression_dir),
        valid_db_genes=valid_db_genes,
        n_specific=top_n_generic,
        n_fallback_scan=top_n_generic,
        filter_mode="protein_coding",
        genome_db=genome_db,
    )

    # --- PHASE 2: Build Sequence Registry ---
    all_target_genes = set().union(*cell_line_top_genes.values(), global_fetch_set)
    click.echo(f"Fetching sequences for {len(all_target_genes)} valid genes...")

    gene_to_data = get_locus_to_data_dict(include_introns=False, gene_subset=list(all_target_genes))

    ref_registry = build_gene_sequence_registry(genes=list(all_target_genes), gene_to_data=gene_to_data)

    # --- PHASE 3: Calculate Weights ---
    cai_weights_map = {}

    # A. Specific Weights
    click.echo("\nCalculating cell-specific weights (Top 300)...")
    for cell_name, genes in cell_line_top_genes.items():
        cds_list = [
            ref_registry[g]["cds"].sequence
            for g in genes
            if g in ref_registry and ref_registry[g].get("cds") and ref_registry[g]["cds"].sequence
        ]
        cds_list = cds_list[:top_n]

        if len(cds_list) == top_n:
            scorer = build_scorer_from_reference(cds_list, pseudocount=pseudocount)
            cai_weights_map[cell_name] = scorer.weights.to_dict()
            click.echo(f"  ✓ {cell_name} weights ready ({len(cds_list)} genes).")
        else:
            click.echo(f"  ⚠ {cell_name} has insufficient sequences.")

    # B. Generic Weights (Consensus Intersection)
    click.echo(f"\nCalculating Generic Weights (Intersection of Top {top_n_generic})...")

    fallback_cds_list = [
        ref_registry[g]["cds"].sequence
        for g in fallback_genes
        if g in ref_registry and ref_registry[g].get("cds") and ref_registry[g]["cds"].sequence
    ]
    fallback_cds_list = fallback_cds_list[:top_n]

    if fallback_cds_list:
        scorer = build_scorer_from_reference(fallback_cds_list)
        cai_weights_map["Generic"] = scorer.weights.to_dict()
        click.echo(f"  ✓ Generic weights ready ({len(fallback_cds_list)} genes).")
    else:
        click.echo("  ⚠ Warning: No sequences found for fallback genes.")

    # --- PHASE 4: Save ---
    with open(out_path, "w") as f:
        json.dump(cai_weights_map, f, indent=4)

    click.echo(click.style(f"\nSUCCESS: CAI weights saved to {out_path}", fg="green"))


@main.command()
@click.option("--force", is_flag=True, help="Force redownload if file exists.")
def setup_mrna_halflife(force):
    """
    Downloads the mRNA half-life dataset (TTDB) from Zenodo and converts it to a
    lean Parquet (only the columns the loader needs) for fast loading.
    """
    # TTDB (Transcriptome Turnover Database) stability dataset — database as of
    # 2026-05-24. Fetched from a Zenodo mirror, verified byte-identical (SHA256
    # ec0c1f90...) to TTDB's official `species_stability_no_threshold.csv.gz`
    # Google Drive download listed at https://sysbio.gzzoc.com/ttdb/download.html.
    # Zenodo is used for a stable, immutable, DOI-backed HTTP download.
    ZENODO_RECORD = "20368324"
    GZ_NAME = "mrna_half_life.csv.gz"
    EXPECTED_SHA256 = "ec0c1f90eed96516de510c484a7a7dd4bbd10d253306710bb33e714f88c5c135"

    data_dir = get_data_dir()
    os.makedirs(data_dir, exist_ok=True)
    gz_path = os.path.join(data_dir, GZ_NAME)
    parquet_path = os.path.join(data_dir, "mrna_half_life.parquet")

    click.echo("Initializing mRNA half-life setup (TTDB via Zenodo)...")
    click.echo(f"Target path: {parquet_path}")

    _ensure_zenodo_table(
        ZENODO_RECORD,
        GZ_NAME,
        gz_path,
        parquet_path,
        EXPECTED_SHA256,
        "sha256",
        force,
        usecols=HALFLIFE_SOURCE_COLUMNS,
    )


@main.command()
@click.option("--force", is_flag=True, help="Force redownload.")
def setup_rrna(force):
    """Download the cytoplasmic rRNA reference FASTA from a pinned Zenodo mirror into
    the data dir, verifying its SHA1.

    Origin: the four human cytoplasmic rRNA RefSeq records -- 18S NR_003286.4, 5.8S NR_003285.3,
    28S NR_003287.4, 5S NR_023363.1 -- fetched once from NCBI nuccore and frozen on Zenodo for
    reproducibility (avoids per-run NCBI rate limits).
    """
    from tauso.features.hybridization.off_target.rrna_targets import REFERENCE_FILENAME, reference_path

    dest = str(reference_path())
    if os.path.exists(dest) and not force and sha1_file(dest) == RRNA_SHA1:
        echo_ok(f"{REFERENCE_FILENAME} exists (SHA1 verified).")
        return
    download_with_progress(
        _zenodo_file_url(ZENODO_RRNA_RECORD, REFERENCE_FILENAME), dest, label=f"    {REFERENCE_FILENAME}"
    )
    verify_hash_or_exit(dest, RRNA_SHA1, algo="sha1")
    echo_ok(f"Downloaded {REFERENCE_FILENAME} (SHA1 verified).")


@main.command(name="setup-model")
@click.option("--version", default=None, help="Model version to fetch (default: the deployed default).")
@click.option("--force", is_flag=True, help="Force redownload if the model file exists.")
def setup_model(version, force):
    """Download the trained ASO-efficacy booster (tauso_score) from Zenodo into
    <data_dir>/models/, verifying its md5, so inference finds it locally instead of
    fetching it on first use. The per-version registry (Zenodo record + md5) lives in
    tauso.inference.predict; this command just provisions it like the other setup-* assets.
    """
    from tauso.inference.predict import DEFAULT_VERSION, ensure_model

    version = version or DEFAULT_VERSION
    click.echo(f"Fetching tauso_score model '{version}' from Zenodo...")
    try:
        dest = ensure_model(version, force=force)
        echo_ok(f"Model ready and verified: {dest}")
    except Exception as e:
        echo_err(f"Error setting up model: {e}")
        sys.exit(1)


@main.command(name="setup-tgcn")
@click.option(
    "--organism",
    type=click.Choice(["human"], case_sensitive=False),
    default="human",
    show_default=True,
    help="Organism whose tGCN to fetch. Available: human (GtRNAdb Hsapi38).",
)
@click.option("--force", is_flag=True, help="Force refetch if file exists.")
def setup_tgcn(organism, force):
    """
    Fetch a tRNA gene copy number (tGCN) table from GtRNAdb and write it
    to TAUSO_DATA_DIR. The tAI feature loads this table at first use;
    without it, populate_tai cannot run.
    """
    from codonbias.scores import fetch_GCN_from_GtRNAdb

    src = TGCNSource[organism.upper()]
    spec = src.value

    data_dir = get_data_dir()
    os.makedirs(data_dir, exist_ok=True)
    destination = os.path.join(data_dir, spec.filename)

    click.echo(f"Initializing {organism} tGCN setup (GtRNAdb {spec.gtrnadb_genome})...")
    click.echo(f"Target path: {destination}")

    if os.path.exists(destination) and not force:
        verify_hash_or_exit(destination, spec.sha256, algo="sha256")
        echo_ok("Existing file matches expected SHA256. Skipping download.")
        return

    try:
        click.echo("Fetching tGCN from GtRNAdb (requires lxml)...")
        df = fetch_GCN_from_GtRNAdb(genome=spec.gtrnadb_genome, domain=spec.gtrnadb_domain)
        df = df.sort_values("anti_codon").reset_index(drop=True)
        df.to_csv(destination, index=False)
        verify_hash_or_exit(destination, spec.sha256, algo="sha256")
        echo_ok(f"Fetched and verified: {destination} ({len(df)} anti-codons).")

    except Exception as e:
        echo_err(f"Error fetching tGCN: {e}")
        sys.exit(1)


@main.command(name="setup-attract")
@click.option("--force", is_flag=True, help="Force redownload if files exist.")
def setup_attract(force):
    """
    Downloads the ATtRACT RBP database files from Zenodo.

    Source: https://zenodo.org/records/20366079 (DOI 10.5281/zenodo.20366079),
    a frozen mirror of the ATtRACT DB. Installs into the 'attract' subfolder of
    TAUSO_DATA_DIR.
    """
    ZENODO_RECORD = "20366079"
    FILES = {
        "RBS_motifs_Homo_sapiens.csv": "12927c0ca4655b4f040ea5e7e3406dc7",
        "pwm.txt": "7ff45b94c14d1b992680f561e0a038a4",
    }

    dest_dir = os.path.join(get_data_dir(), "attract")
    os.makedirs(dest_dir, exist_ok=True)
    click.echo(f"Target directory: {dest_dir}")

    for name, expected_md5 in FILES.items():
        destination = os.path.join(dest_dir, name)
        _ensure_zenodo_content_file(ZENODO_RECORD, name, destination, expected_md5, "md5", force)


_RIBOSEQ_ZENODO_RECORD = "20435808"
_RIBOSEQ_TRACKS = (
    # (filename, expected_md5)
    ("human_unselected_40S.RiboProElong.bw", "c1a06bf87fbee3d66f8422922dddd709"),
    ("human_unselected_80S.RiboCov.bw", "ca2ef00c545254bc0c3eceecc58fd2a2"),
)


@main.command(name="setup-riboseq")
@click.option("--force", is_flag=True, help="Force redownload if file exists.")
def setup_riboseq(force):
    """Downloads the Wagner 2020 ribo-seq bigWigs into TAUSO_DATA_DIR.

    40S unselected scanning (~7 MB) + 80S unselected elongation (~13 MB), both
    from HEK293T Sel-TCP-seq (GEO GSE139131) and mirrored at Zenodo record
    10.5281/zenodo.20435808. The features built from these are gene-level
    translation proxies, not cell-context: the same value applies to every
    TAUSO row regardless of cell line.
    """
    data_dir = get_data_dir()
    os.makedirs(data_dir, exist_ok=True)

    for filename, expected_md5 in _RIBOSEQ_TRACKS:
        destination = os.path.join(data_dir, filename)
        click.echo(f"Target path: {destination}")
        _ensure_zenodo_content_file(_RIBOSEQ_ZENODO_RECORD, filename, destination, expected_md5, "md5", force)


@main.command(name="setup-features")
@click.option("--run", default="oligo", show_default=True, help="Feature-set run to download.")
@click.option("--force", is_flag=True, help="Force redownload even if present.")
@click.option(
    "--include-competition",
    is_flag=True,
    help="Also download the competitor-score shards (OligoWalk/miRanda/sFold/PFRED/OligoAI).",
)
def setup_features(run, force, include_competition):
    """
    Download the frozen feature cache (one wide parquet) from Zenodo into
    TAUSO_DATA_DIR/features/<run>/.

    This is the index-keyed cache for the original ASO dataset -- an OPTIONAL fast path so
    modeling/notebooks don't recompute the full feature matrix. Computing features for new
    inputs goes through the calculators and does not need it. The cache sits in the store
    folder alongside any per-feature dev shards.

    With --include-competition, also fetches the competitor-score shards into
    features/<run>/_patches/ -- the opt-in baselines loaded only when load_competition=True.
    """
    from tauso.populate.feature_cache import ensure_cache, ensure_competition_cache

    try:
        dest = ensure_cache(run, force=force)
    except FileNotFoundError as e:
        echo_err(str(e))
        sys.exit(1)
    echo_ok(f"Feature cache ready: {dest}")

    if include_competition:
        paths = ensure_competition_cache(run, force=force)
        echo_ok(f"Competition shards ready: {len(paths)} files")


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

            gffutils.create_db(
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

        # Pre-build the coordinate-only locus pickle so populate runs (and SLURM
        # tasks) skip the ~24s gff traversal and load it from cache instead.
        click.echo("Building locus coordinate cache...")
        try:
            cache_path = build_locus_cache(genome=genome, include_introns=True, canonical_only=True)
            click.echo(click.style(f"✓ Locus cache built at {cache_path}", fg="green"))
        except Exception as e:
            # Non-fatal: the cache is an optimization; populate falls back to the gff db.
            click.echo(click.style(f"⚠ Could not build locus cache ({e}); will build on first use.", fg="yellow"))


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
    sequence = rna_to_dna(sequence)

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
@click.option(
    "--march",
    default=None,
    help="-march for the raccess build (default: native, or $RACCESS_MARCH). "
    "Use 'x86-64-v2' for a portable binary that runs across heterogeneous "
    "cluster CPUs (native can crash with SIGILL when build and run nodes differ).",
)
def setup_raccess(ctx, force_clone, march):
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

    # --march flag wins, else inherit $RACCESS_MARCH, else the script's 'native' default
    env = os.environ.copy()
    if march:
        env["RACCESS_MARCH"] = march

    click.echo(f"+ {' '.join(cmd)} (RACCESS_MARCH={env.get('RACCESS_MARCH', 'native')})")
    try:
        subprocess.run(cmd, check=True, env=env)
    except subprocess.CalledProcessError:
        logger.error("install_raccess.sh failed, consider installing zlib1g-dev")
        sys.exit(1)


if __name__ == "__main__":
    main()
