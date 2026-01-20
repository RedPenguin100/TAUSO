import os
import sys
import shutil
import logging
import subprocess
import click
import glob
import time
import itertools
import pandas as pd
from tauso.data.data import get_paths, load_db

# Setup logger
logger = logging.getLogger(__name__)
if not logger.handlers:
    handler = logging.StreamHandler(sys.stderr)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)


def get_bowtie_index_base(genome="GRCh38", force_rebuild=False):
    """
    Returns the base path for the Bowtie index of the specific genome.
    """
    if not shutil.which("bowtie"):
        raise RuntimeError("Bowtie binary not found. Please install: 'micromamba install -c bioconda bowtie'")

    paths = get_paths(genome)
    fasta_path = paths['fasta']

    # Check if FASTA exists first
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Genome FASTA not found: {fasta_path}. Please run setup-genome.")

    # Unique index folder per genome
    index_dir = os.path.join(os.path.dirname(fasta_path), f"{genome}_bowtie_index")
    os.makedirs(index_dir, exist_ok=True)

    index_base = os.path.join(index_dir, "genome")
    sentinel_file = os.path.join(index_dir, "SUCCESS")

    std_exts = [".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt"]
    lrg_exts = [".1.ebwtl", ".2.ebwtl", ".3.ebwtl", ".4.ebwtl", ".rev.1.ebwtl", ".rev.2.ebwtl"]

    def index_status(base):
        if not os.path.exists(sentinel_file):
            if glob.glob(f"{base}*.ebwt*"): return 'partial'
            return 'missing'

        std_exists = [os.path.exists(base + ext) for ext in std_exts]
        lrg_exists = [os.path.exists(base + ext) for ext in lrg_exts]

        if all(std_exists) or all(lrg_exists): return 'complete'
        return 'partial'

    def clean_index(base):
        for f in glob.glob(f"{base}*.ebwt*"):
            try:
                os.remove(f)
            except OSError:
                pass
        if os.path.exists(sentinel_file):
            try:
                os.remove(sentinel_file)
            except OSError:
                pass

    status = index_status(index_base)

    if not force_rebuild:
        if status == 'complete':
            logger.debug(f"Found valid Bowtie index: {index_base}")
            return index_base
        elif status == 'partial':
            logger.warning(f"Found partial/corrupt index for {genome}. Rebuilding...")
            clean_index(index_base)
    elif force_rebuild and status != 'missing':
        logger.info(f"Force rebuild requested for {genome}...")
        clean_index(index_base)

    logger.info(f"Building Bowtie index for {genome}...")
    try:
        cmd = ["bowtie-build", fasta_path, index_base]

        with subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True) as proc:
            spinner = itertools.cycle(['-', '/', '|', '\\'])
            while proc.poll() is None:
                sys.stdout.write(f"\r  Compiling {genome} Index... {next(spinner)}")
                sys.stdout.flush()
                time.sleep(0.5)

            sys.stdout.write(f"\r  Compiling {genome} Index... Done!      \n")

            if proc.returncode != 0:
                stderr_output = proc.stderr.read()
                raise subprocess.CalledProcessError(proc.returncode, cmd, stderr_output)

        with open(sentinel_file, "w") as f:
            f.write("Index build successful.")
        logger.info(f"Index built: {index_base}")

    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to build index: {e.stderr}")
        clean_index(index_base)
        raise RuntimeError("Bowtie indexing failed.")
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        clean_index(index_base)
        raise e

    return index_base


def run_bowtie_search(sequence, genome="GRCh38", max_mismatches=3):
    """
    Runs Bowtie 1 alignment.
    Returns:
        hits_list: List of all hit dictionaries (for annotation)
        counts_dict: Dictionary of counts {'mismatches0': X, 'mismatches1': Y...}
    """
    index_base = get_bowtie_index_base(genome=genome)

    cmd = [
        "bowtie",
        "-v", str(max_mismatches),
        "-a",  # Report all valid alignments
        "-S", "--sam-nohead",
        index_base,
        "-c", sequence
    ]

    try:
        process = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Bowtie search failed: {e.stderr}")
        return [], {f'mismatches{i}': 0 for i in range(max_mismatches + 1)}

    hits = []

    # --- IMPLEMENTATION OF YOUR REQUEST ---
    # Initialize the columns/counters you wanted
    counts = {f'mismatches{i}': 0 for i in range(max_mismatches + 1)}

    for line in process.stdout.splitlines():
        if not line.strip(): continue
        parts = line.split('\t')
        flag = int(parts[1])
        if flag & 4: continue  # Unmapped

        chrom = parts[2]
        start_pos = int(parts[3]) - 1

        # Parse mismatches from NM tag
        mismatches = 0
        for tag in parts[11:]:
            if tag.startswith("NM:i:"):
                mismatches = int(tag.split(":")[2])
                break

        # 1. Do the ++ for the specific mismatch column
        if mismatches <= max_mismatches:
            counts[f'mismatches{mismatches}'] += 1

        # 2. Keep the hit data (needed for 'annotate_hits' later)
        hits.append({
            'chrom': chrom,
            'start': start_pos,
            'end': start_pos + len(sequence),
            'strand': '-' if (flag & 16) else '+',
            'mismatches': mismatches,
            'sequence': sequence
        })

    return hits, counts


def annotate_hits(hits_list, genome="GRCh38"):
    """
    Annotates hits using the specific genome database.
    """
    if not hits_list:
        return pd.DataFrame()

    db = load_db(genome=genome)
    annotated = []

    for hit in hits_list:
        chrom = hit['chrom']
        start = hit['start']
        end = hit['end']

        try:
            features = list(db.region(region=(chrom, start, end)))
        except Exception:
            features = []

        gene_id = None
        gene_name = None
        feature_type = "Intergenic"
        current_priority = 0

        for feat in features:
            f_type = feat.featuretype

            # --- UPDATED PRIORITY LOGIC ---
            if f_type == 'exon':
                priority = 4
            elif f_type == 'CDS':  # Added for Bacteria/Yeast compatibility
                priority = 4
            elif f_type == 'intron':
                priority = 2
            elif f_type == 'gene':
                priority = 1
            else:
                priority = 0
            # ------------------------------

            if priority > current_priority:
                current_priority = priority
                feature_type = f_type
                # Ensembl sometimes uses 'gene_name', sometimes 'Name', sometimes just 'gene_id'
                # This fallback chain covers Human (gene_name) and Bacteria (Name)
                gene_name = feat.attributes.get('gene_name',
                                                feat.attributes.get('Name',
                                                                    feat.attributes.get('gene_id', [None])))[0]

                gene_id = feat.attributes.get('gene_id', [None])[0]

        hit_copy = hit.copy()
        hit_copy['gene_id'] = gene_id
        hit_copy['gene_name'] = gene_name
        hit_copy['region_type'] = feature_type
        annotated.append(hit_copy)

    return pd.DataFrame(annotated)


def find_all_gene_off_targets(sequence, genome="GRCh38", max_mismatches=3):
    """
    Main entry point for CLI.
    """
    # Unpack the tuple (hits, counts)
    hits_list, counts_dict = run_bowtie_search(sequence, genome=genome, max_mismatches=max_mismatches)

    # Use the hits list for annotation as before
    df = annotate_hits(hits_list, genome=genome)

    # Optional: You could print the counts_dict here if debugging
    # print(f"Counts: {counts_dict}")

    return df