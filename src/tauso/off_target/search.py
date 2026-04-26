import glob
import itertools
import logging
import os
import shutil
import subprocess
import sys
import time

import pandas as pd
import pyranges as pr

from ..data.data import get_paths, load_db, load_pyranges_gtf
from ..timer import Timer

# TODO: standardize logger in the entire package
# Setup logger
logger = logging.getLogger(__name__)
if not logger.handlers:
    handler = logging.StreamHandler(sys.stderr)
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)


def get_bowtie_index_base(genome="GRCh38", force_rebuild=False, threads=1, mem_per_thread_mb=800):
    """
    Returns the base path for the Bowtie index of the specific genome.
    """
    if not shutil.which("bowtie"):
        raise RuntimeError("Bowtie binary not found. Please install: 'micromamba install -c bioconda bowtie'")

    paths = get_paths(genome)
    fasta_path = paths["fasta"]

    # Check if FASTA exists first
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Genome FASTA not found: {fasta_path}. Please run setup-genome.")

    # Unique index folder per genome
    index_dir = os.path.join(os.path.dirname(fasta_path), f"{genome}_bowtie_index")
    os.makedirs(index_dir, exist_ok=True)

    index_base = os.path.join(index_dir, "genome")
    sentinel_file = os.path.join(index_dir, "SUCCESS")

    std_exts = [
        ".1.ebwt",
        ".2.ebwt",
        ".3.ebwt",
        ".4.ebwt",
        ".rev.1.ebwt",
        ".rev.2.ebwt",
    ]
    lrg_exts = [
        ".1.ebwtl",
        ".2.ebwtl",
        ".3.ebwtl",
        ".4.ebwtl",
        ".rev.1.ebwtl",
        ".rev.2.ebwtl",
    ]

    def index_status(base):
        if not os.path.exists(sentinel_file):
            if glob.glob(f"{base}*.ebwt*"):
                return "partial"
            return "missing"

        std_exists = [os.path.exists(base + ext) for ext in std_exts]
        lrg_exists = [os.path.exists(base + ext) for ext in lrg_exts]

        if all(std_exists) or all(lrg_exists):
            return "complete"
        return "partial"

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
        if status == "complete":
            logger.debug(f"Found valid Bowtie index: {index_base}")
            return index_base
        elif status == "partial":
            logger.warning(f"Found partial/corrupt index for {genome}. Rebuilding...")
            clean_index(index_base)
    elif force_rebuild and status != "missing":
        logger.info(f"Force rebuild requested for {genome}...")
        clean_index(index_base)

    logger.info(f"Building Bowtie index for {genome}...")
    log_file_path = os.path.join(index_dir, "bowtie_build.log")
    logger.info(f"Detailed build logs will be saved to: {log_file_path}")

    try:
        # 1 MB = 1,048,576 bytes.
        # In Bowtie's BWT construction, 1 suffix pointer requires ~4 bytes of RAM.
        max_bytes = mem_per_thread_mb * 1024 * 1024
        bmax_suffixes = int(max_bytes / 4)

        cmd = [
            "bowtie-build",
            "--threads",
            str(threads),
            "--packed",
            "--bmax",
            str(bmax_suffixes),
            "--dcv",
            "2048",
            "--offrate",
            "6",
            fasta_path,
            index_base,
        ]

        # Route both stdout and stderr to the log file instead of DEVNULL
        with open(log_file_path, "w") as log_file:
            with subprocess.Popen(cmd, stdout=log_file, stderr=subprocess.STDOUT, text=True) as proc:
                spinner = itertools.cycle(["-", "/", "|", "\\"])
                while proc.poll() is None:
                    sys.stdout.write(f"\r  Compiling {genome} Index... {next(spinner)}")
                    sys.stdout.flush()
                    time.sleep(0.5)

            sys.stdout.write(f"\r  Compiling {genome} Index... Done!      \n")

            if proc.returncode != 0:
                # Grab the last 15 lines of the log to see the actual error
                with open(log_file_path, "r") as f:
                    log_lines = f.readlines()
                    error_tail = "".join(log_lines[-15:]) if log_lines else "Log empty. Process killed by OS?"

                raise subprocess.CalledProcessError(proc.returncode, cmd, error_tail)

        with open(sentinel_file, "w") as f:
            f.write("Index build successful.")
        logger.info(f"Index built: {index_base}")

    except subprocess.CalledProcessError as e:
        logger.error(
            f"Failed to build index. Exit Code: {e.returncode}\n--- Last Log Output ---\n{e.stderr}\n-----------------------"
        )
        clean_index(index_base)
        raise RuntimeError(f"Bowtie indexing failed. See {log_file_path} for details.")
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
        "-v",
        str(max_mismatches),
        "-a",  # Report all valid alignments
        "-S",
        "--sam-nohead",
        "-x",  # <--- ADD THIS LINE
        index_base,
        "-c",
        sequence,
    ]

    try:
        process = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Bowtie search failed: {e.stderr}")
        return [], {f"mismatches{i}": 0 for i in range(max_mismatches + 1)}

    hits = []

    # --- IMPLEMENTATION OF YOUR REQUEST ---
    # Initialize the columns/counters you wanted
    counts = {f"mismatches{i}": 0 for i in range(max_mismatches + 1)}

    for line in process.stdout.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        flag = int(parts[1])
        if flag & 4:
            continue  # Unmapped

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
            counts[f"mismatches{mismatches}"] += 1

        # 2. Keep the hit data (needed for 'annotate_hits' later)
        hits.append(
            {
                "chrom": chrom,
                "start": start_pos,
                "end": start_pos + len(sequence),
                "strand": "-" if (flag & 16) else "+",
                "mismatches": mismatches,
                "sequence": sequence,
            }
        )

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
        chrom = hit["chrom"]
        start = hit["start"]
        end = hit["end"]

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
            if f_type == "exon":
                priority = 4
            elif f_type == "CDS":  # Added for Bacteria/Yeast compatibility
                priority = 4
            elif f_type == "intron":
                priority = 2
            elif f_type == "gene":
                priority = 1
            else:
                priority = 0
            # ------------------------------

            if priority > current_priority:
                current_priority = priority
                feature_type = f_type
                # Ensembl sometimes uses 'gene_name', sometimes 'Name', sometimes just 'gene_id'
                # This fallback chain covers Human (gene_name) and Bacteria (Name)
                gene_name = feat.attributes.get(
                    "gene_name",
                    feat.attributes.get("Name", feat.attributes.get("gene_id", [None])),
                )[0]

                gene_id = feat.attributes.get("gene_id", [None])[0]

        hit_copy = hit.copy()
        hit_copy["gene_id"] = gene_id
        hit_copy["gene_name"] = gene_name
        hit_copy["region_type"] = feature_type
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

    # print(f"Counts: {counts_dict}")
    # print(f"hits_list: {hits_list}")

    return df


def run_bowtie_search_bulk(fasta_path, genome="GRCh38", max_mismatches=0, threads=16):
    """
    Runs Bowtie 1 alignment on a whole FASTA file.
    """
    index_base = get_bowtie_index_base(genome=genome)
    sam_output = fasta_path.replace(".fasta", ".sam")

    cmd = [
        "bowtie",
        "-v",
        str(max_mismatches),
        "-a",  # Report all valid alignments
        "-S",
        "--sam-nohead",
        "-p",
        str(threads),  # MULTITHREADING ENABLED
        "-f",  # FASTA INPUT (changed from -c)
        "-x",
        index_base,
        fasta_path,  # Input file
        sam_output,  # Output file
    ]

    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Bowtie bulk search failed: {e.stderr}")
        return []

    hits = []

    # Parse the massive SAM output line by line
    with open(sam_output, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("@"):
                continue

            parts = line.split("\t")
            flag = int(parts[1])
            if flag & 4:
                continue  # Unmapped

            # Because we wrote the FASTA with >SEQUENCE, parts[0] is the sequence itself
            sequence_id = parts[0]
            chrom = parts[2]
            start_pos = int(parts[3]) - 1

            hits.append(
                {"sequence": sequence_id, "chrom": chrom, "start": start_pos, "end": start_pos + len(sequence_id)}
            )

    if os.path.exists(sam_output):
        os.remove(sam_output)

    return hits


def annotate_hits_bulk(hits_list, genome):
    """
    Lightning-fast pyranges intersection replacing the SQLite loop.
    Pure spatial intersection: returns ALL genes an ASO touches,
    ignoring feature types, introns, or biotypes.
    """
    if not hits_list:
        return {}

    # Assumes get_paths is defined globally in your script
    gr_genome = load_pyranges_gtf(get_paths(genome)["gtf"])

    # 1. Convert your raw Bowtie hits to a Pandas DataFrame
    df_hits = pd.DataFrame(hits_list)

    # Pyranges explicitly requires these exact column names (capitalized)
    df_hits = df_hits.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})

    # Shrink memory footprint drastically
    df_hits["Chromosome"] = df_hits["Chromosome"].astype("category")
    df_hits["sequence"] = df_hits["sequence"].astype("category")  # If sequence strings are repetitive
    df_hits["Start"] = pd.to_numeric(df_hits["Start"], downcast="integer")
    df_hits["End"] = pd.to_numeric(df_hits["End"], downcast="integer")

    # Convert to a PyRanges object
    # gr_hits = pr.PyRanges(df_hits)

    with Timer(f"Intersecting {len(df_hits)} hits against the genome in chunks..."):
        chunk_size = 1_000_000
        df_res_list = []

        for i in range(0, len(df_hits), chunk_size):
            print(f"  Processing chunk {i} to {i + chunk_size}...")

            # 1. Slice the dataframe
            df_chunk = df_hits.iloc[i : i + chunk_size]

            # 2. Convert chunk to PyRanges
            gr_chunk = pr.PyRanges(df_chunk)

            # 3. Intersect just this chunk
            intersected_chunk = gr_chunk.join(gr_genome, apply_strand_suffix=False)

            # 4. Save the resulting dataframe and free memory
            if not intersected_chunk.df.empty:
                df_res_list.append(intersected_chunk.df)

    if not df_res_list:
        return {}

    # Combine all the chunked results back into one dataframe
    df_res = pd.concat(df_res_list, ignore_index=True)

    # df_res = intersected.df

    if df_res.empty:
        return {}

    # 3. Resolve the gene name based on the fallback chain
    if "gene_name" in df_res.columns:
        df_res["resolved_gene_name"] = df_res["gene_name"]
    elif "Name" in df_res.columns:
        df_res["resolved_gene_name"] = df_res["Name"]
    elif "gene_id" in df_res.columns:
        df_res["resolved_gene_name"] = df_res["gene_id"]
    else:
        df_res["resolved_gene_name"] = None

    # Drop anything that didn't map to a name
    df_res.dropna(subset=["resolved_gene_name"], inplace=True)

    # 4. AGGREGATE BACK TO DICTIONARY (No sorting or dropping duplicates needed)
    # Group by the ASO sequence and collect all unique gene names it hit
    seq_to_genes_series = df_res.groupby("sequence", observed=True)["resolved_gene_name"].unique()
    # Convert pandas Series of arrays to a standard Python dictionary of lists
    seq_to_genes = {seq: list(genes) for seq, genes in seq_to_genes_series.items()}

    return seq_to_genes


def find_all_gene_off_targets_BULK(fasta_path, genome="GRCh38", threads=16, max_mismatches=0):
    """
    Main bulk entry point. Replaces your old find_all_gene_off_targets.
    """
    with Timer("Running bowtie"):
        hits_list = run_bowtie_search_bulk(fasta_path, genome=genome, max_mismatches=max_mismatches, threads=threads)

    with Timer("Annotate bowtie hits"):
        # 2. Annotate the hits in bulk and get the mapping dictionary
        seq_to_genes = annotate_hits_bulk(hits_list, genome=genome)

    return seq_to_genes
