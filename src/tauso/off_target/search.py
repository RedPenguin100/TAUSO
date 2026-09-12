import glob
import logging
import os
import shutil
import subprocess
import tempfile
import time
from functools import lru_cache

import numpy as np
import pandas as pd
import pyranges as pr
from ncls import NCLS

from ..data.data import ANNOTATION_PRIORITY, get_paths, load_gene_intervals, load_gtf_pyranges_gene_only
from ..debug import log_memory_usage
from ..timer import Timer

logger = logging.getLogger(__name__)


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
                start_time = time.time()
                while proc.poll() is None:
                    elapsed = time.time() - start_time
                    index_files = glob.glob(f"{index_base}*.ebwt*")
                    total_gb = sum(os.path.getsize(f) for f in index_files if os.path.exists(f)) / (1024**3)
                    logger.info(
                        f"Building {genome} index... {total_gb:.2f} GB written | elapsed: {int(elapsed // 60)}m{int(elapsed % 60)}s"
                    )
                    time.sleep(30)

            logger.info(f"Compiling {genome} Index... Done!")

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
    Raises RuntimeError on alignment failure, including Bowtie's diagnostics.
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
        "-x",
        index_base,
        "-c",
        sequence,
    ]

    try:
        process = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Bowtie search failed: {e}\n{e.stderr or ''}") from e

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


@lru_cache(maxsize=2)
def _annotation_index(genome):
    """The ranked annotation as (column arrays, one interval tree per chrom/strand).

    One tree per (chrom, strand) so a hit only ever searches the strand it can be antisense
    to. NCLS is half-open and the GTF is closed, so ends are widened by one; querying with
    `end + 1` as well then reproduces the closed-interval overlap gffutils tests for.
    """
    df = load_gene_intervals(genome)
    columns = {
        "start": df["start"].to_numpy(np.int64),
        "end": df["end"].to_numpy(np.int64),
        "priority": df["featuretype"].map(ANNOTATION_PRIORITY).to_numpy(np.int16),
        "featuretype": df["featuretype"].to_numpy(object),
        "gene_name": df["gene_name"].to_numpy(object),
        "gene_id": df["gene_id"].to_numpy(object),
    }
    trees = {}
    for key, sub in df.groupby(["chrom", "strand"], sort=False, observed=True):
        rows = sub.index.to_numpy(np.int64)
        trees[key] = NCLS(columns["start"][rows], columns["end"][rows] + 1, rows)
    return columns, trees


def annotate_hits(hits_list, genome="GRCh38"):
    """Annotate each hit with its gene and region, counting a gene only when the ASO is antisense to
    it (hit strand opposite the gene's). Same-strand overlaps are ignored, so the hit is Intergenic.
    Annotation loading and interval query errors propagate to the caller.

    Every hit is resolved in one pass over an in-memory interval index. Querying the annotation
    per hit instead costs a full chromosome scan each time, because the overlap test cannot use
    the database's index -- minutes rather than seconds once a low-complexity ASO aligns widely.
    """
    if not hits_list:
        return pd.DataFrame()

    columns, trees = _annotation_index(genome)
    hits = pd.DataFrame(hits_list)
    chroms = hits["chrom"].to_numpy(object)
    starts = hits["start"].to_numpy(np.int64)
    ends = hits["end"].to_numpy(np.int64)
    # A hit only annotates against genes it is antisense to, so search the opposite strand.
    antisense = np.where(hits["strand"].to_numpy(object) == "+", "-", "+")

    hit_rows, feature_rows = [], []
    for key, group in pd.Series(np.arange(len(hits))).groupby([chroms, antisense], sort=False):
        tree = trees.get(key)
        if tree is None:
            continue  # a contig the annotation does not cover
        rows = group.to_numpy()
        local, features = tree.all_overlaps_both(starts[rows], ends[rows] + 1, np.arange(len(rows), dtype=np.int64))
        if len(local):
            hit_rows.append(rows[local])
            feature_rows.append(features)

    gene_id = np.full(len(hits), None, dtype=object)
    gene_name = np.full(len(hits), None, dtype=object)
    region_type = np.full(len(hits), "Intergenic", dtype=object)

    if hit_rows:
        hit_row = np.concatenate(hit_rows)
        feature_row = np.concatenate(feature_rows)
        # Highest priority wins; ties go to the earliest feature, ordered as the annotation is.
        order = np.lexsort(
            (
                feature_row,
                columns["end"][feature_row],
                columns["start"][feature_row],
                -columns["priority"][feature_row],
                hit_row,
            )
        )
        hit_row, feature_row = hit_row[order], feature_row[order]
        best = np.flatnonzero(np.r_[True, hit_row[1:] != hit_row[:-1]])
        winner_hit, winner_feature = hit_row[best], feature_row[best]
        region_type[winner_hit] = columns["featuretype"][winner_feature]
        gene_name[winner_hit] = columns["gene_name"][winner_feature]
        gene_id[winner_hit] = columns["gene_id"][winner_feature]

    hits["gene_id"] = gene_id
    hits["gene_name"] = gene_name
    hits["region_type"] = region_type
    return hits


def find_all_gene_off_targets(sequence, genome="GRCh38", max_mismatches=3):
    """
    Main entry point for CLI.
    """
    # Unpack the tuple (hits, counts)
    hits_list, counts_dict = run_bowtie_search(sequence, genome=genome, max_mismatches=max_mismatches)

    # Use the hits list for annotation as before
    df = annotate_hits(hits_list, genome=genome)

    # logger.debug(f"Counts: {counts_dict}")
    # logger.debug(f"hits_list: {hits_list}")

    return df


def run_bowtie_search_bulk(fasta_path, genome="GRCh38", max_mismatches=0, threads=16):
    """
    Runs Bowtie 1 alignment on a whole FASTA file.
    Raises RuntimeError on alignment failure, including Bowtie's diagnostics.
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
        raise RuntimeError(f"Bowtie bulk search failed: {e}\n{e.stderr or ''}") from e

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

            # Because we wrote the FASTA with >ASO_SEQUENCE, parts[0] is the sequence itself
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

    gr_genome = load_gtf_pyranges_gene_only(get_paths(genome)["gtf_gz"])

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
            logger.debug("Processing chunk %d to %d...", i, i + chunk_size)

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


@log_memory_usage
def find_all_gene_off_targets_BULK(fasta_path, genome="GRCh38", threads=16, max_mismatches=0):
    logger.debug("[Find_OT] Running bowtie")
    hits_list = run_bowtie_search_bulk(fasta_path, genome=genome, max_mismatches=max_mismatches, threads=threads)

    logger.debug("[Find_OT] Annotate hits")
    # 2. Annotate the hits in bulk and get the mapping dictionary
    seq_to_genes = annotate_hits_bulk(hits_list, genome=genome)

    return seq_to_genes


def find_all_gene_off_targets_bulk_sequences(sequences, genome="GRCh38", threads=16, max_mismatches=0):
    """Genes each sequence aligns to, in a single Bowtie pass. Returns ``{sequence: [gene, ...]}``.

    Sequence-taking counterpart to `find_all_gene_off_targets_BULK`; owns the temp FASTA so callers
    holding sequences rather than a file do not have to.
    """
    seqs = list(dict.fromkeys(sequences))
    if not seqs:
        return {}

    with tempfile.TemporaryDirectory() as work:
        fasta_path = os.path.join(work, "sequences.fasta")
        with open(fasta_path, "w") as f:
            for s in seqs:
                f.write(f">{s}\n{s}\n")
        return find_all_gene_off_targets_BULK(fasta_path, genome=genome, threads=threads, max_mismatches=max_mismatches)


def count_offtarget_matches_bulk(sequences, genome="GRCh38", max_mismatches=2, threads=16, exclude_regions=None):
    """Count genome matches per sequence at each mismatch distance, in a single Bowtie pass.

    Aligns every sequence to `genome` (both strands, all alignments up to `max_mismatches`) and tallies,
    for each sequence, how many genomic loci it matches at exactly k mismatches (k in 0..max_mismatches).
    Returns ``{sequence: {0: n0, 1: n1, ...}}``.
    Raises RuntimeError on alignment failure, including Bowtie's diagnostics.

    `exclude_regions` is an optional iterable of ``(chrom, start, end)`` genomic intervals (0-based,
    half-open); any hit overlapping one is not counted. Pass the on-target gene's locus so the intended
    site and any intragenic near-matches do not inflate the off-target tallies.

    Assumes the Bowtie index exists (see `get_bowtie_index_base`); callers that must not trigger a
    multi-GB index build should check the index SUCCESS sentinel before calling.
    """
    index_base = get_bowtie_index_base(genome=genome)

    seqs = list(dict.fromkeys(sequences))
    counts = {s: {i: 0 for i in range(max_mismatches + 1)} for s in seqs}
    if not seqs:
        return counts

    exclude = list(exclude_regions or [])

    with tempfile.TemporaryDirectory() as work:
        fasta_path = os.path.join(work, "aso.fasta")
        sam_path = os.path.join(work, "aso.sam")
        with open(fasta_path, "w") as f:
            for s in seqs:
                f.write(f">{s}\n{s}\n")

        cmd = [
            "bowtie",
            "-v",
            str(max_mismatches),
            "-a",  # Report all valid alignments
            "-S",
            "--sam-nohead",
            "-p",
            str(threads),
            "-f",  # FASTA input; read name == the sequence
            "-x",
            index_base,
            fasta_path,
            sam_path,
        ]
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Bowtie bulk count failed: {e}\n{e.stderr or ''}") from e

        with open(sam_path) as sam:
            for line in sam:
                if not line.strip() or line.startswith("@"):
                    continue
                parts = line.split("\t")
                flag = int(parts[1])
                if flag & 4:
                    continue  # Unmapped

                seq_id = parts[0]
                if seq_id not in counts:
                    continue

                if exclude:
                    start_pos = int(parts[3]) - 1
                    end_pos = start_pos + len(seq_id)
                    if any(parts[2] == ec and start_pos < ee and es < end_pos for ec, es, ee in exclude):
                        continue

                mismatches = 0
                for tag in parts[11:]:
                    if tag.startswith("NM:i:"):
                        mismatches = int(tag.split(":")[2])
                        break

                if mismatches <= max_mismatches:
                    counts[seq_id][mismatches] += 1

    return counts
