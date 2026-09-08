import gzip
import logging
import os
import tempfile
from functools import lru_cache

import gffutils
import pandas as pd
import pyranges as pr
from platformdirs import user_data_dir
from pyfaidx import Fasta

APP_NAME = "tauso"


logger = logging.getLogger(__name__)


def get_data_dir():
    # 1. Check for Docker/Render Environment Variable override
    env_path = os.environ.get("TAUSO_DATA_DIR")

    if env_path:
        d = env_path
    else:
        # 2. Fallback to standard User Data Dir (e.g., ~/.local/share/tauso)
        d = user_data_dir(APP_NAME)

    # Ensure the directory exists (works for both paths)
    os.makedirs(d, exist_ok=True)
    return d


def get_paths(genome="GRCh38"):
    """
    Returns paths for a specific genome assembly.
    Default: GRCh38
    """
    d = get_data_dir()
    return {
        "dir": d,
        "fasta": os.path.join(d, f"{genome}.fa"),
        "gtf": os.path.join(d, f"{genome}.gtf"),
        "gtf_gz": os.path.join(d, f"{genome}.gtf.gz"),
        "gff": os.path.join(d, f"{genome}.gff3"),
        "gff_gz": os.path.join(d, f"{genome}.gff3.gz"),
        "gff_db": os.path.join(d, f"{genome}.gff3.db"),
        "gtf_db": os.path.join(d, f"{genome}.gtf.db"),
        "gene_intervals": os.path.join(d, f"{genome}.intervals.parquet"),
    }


# Feature types `annotate_hits` ranks a hit against, and the priority it gives each. Introns
# are listed for annotations that carry them; GTFs generally do not (gffutils does not infer
# them either), in which case an intronic hit falls through to its gene.
ANNOTATION_PRIORITY = {"exon": 4, "CDS": 4, "intron": 2, "gene": 1}


def _parse_gene_intervals(gtf_gz_path: str):
    """Read the annotation's ranked features out of the GTF into a flat interval table."""
    df = pd.read_csv(
        gtf_gz_path,
        sep="\t",
        comment="#",
        header=None,
        usecols=[0, 2, 3, 4, 6, 8],
        names=["chrom", "featuretype", "start", "end", "strand", "attributes"],
    )
    # Row order stays the GTF's own line order, which is the order gffutils inserted with and
    # the order `annotate_hits` breaks priority ties on, so it has to survive the round trip.
    df = df[df["featuretype"].isin(ANNOTATION_PRIORITY)].reset_index(drop=True)
    df["gene_id"] = df["attributes"].str.extract(r'gene_id "([^"]+)"', expand=False)
    # Ensembl writes 'gene_name', some annotations only 'Name', and a few neither: fall back
    # down that chain to gene_id, as the per-feature attribute lookup this replaced did.
    name = df["attributes"].str.extract(r'gene_name "([^"]+)"', expand=False)
    name = name.fillna(df["attributes"].str.extract(r'(?:^|; )Name "([^"]+)"', expand=False))
    df["gene_name"] = name.fillna(df["gene_id"])
    df = df.drop(columns="attributes")
    # Missing values reach callers as None, not NaN, matching the attribute lookup's default.
    for column in ("gene_name", "gene_id"):
        df[column] = df[column].astype(object).where(df[column].notna(), None)
    return df


@lru_cache(maxsize=2)
def load_gene_intervals(genome="GRCh38"):
    """Ranked annotation features as a flat interval table, one row per exon/CDS/gene.

    Built from the GTF on first use and cached beside it as parquet (~16 MB, a few seconds),
    then reloaded in well under a second. Cached per process: callers get a shared frame and
    must not mutate it.
    """
    paths = get_paths(genome)
    cache = paths["gene_intervals"]
    if os.path.exists(cache):
        return pd.read_parquet(cache)

    if not os.path.exists(paths["gtf_gz"]):
        raise FileNotFoundError(f"GTF for {genome} not found. Run 'tauso setup-genome --genome {genome}'")
    logger.info("Building the gene interval cache for %s (first use only)...", genome)
    df = _parse_gene_intervals(paths["gtf_gz"])
    tmp = f"{cache}.{os.getpid()}.tmp"  # write-then-rename: a crash must not leave a half file
    df.to_parquet(tmp, index=False)
    os.replace(tmp, cache)
    logger.info("Cached %d annotation features to %s", len(df), cache)
    return df


def load_gtf_db(genome="GRCh38"):
    """Returns the gffutils database for the specified genome."""
    paths = get_paths(genome)
    if not os.path.exists(paths["gtf_db"]):
        raise FileNotFoundError(f"Database for {genome} not found. Run 'tauso setup-genome --genome {genome}'")
    return gffutils.FeatureDB(paths["gtf_db"])


def load_gff_db(genome="GRCh38"):
    """Returns the gffutils database for the specified genome."""
    paths = get_paths(genome)
    if not os.path.exists(paths["gff_db"]):
        raise FileNotFoundError(f"Database for {genome} not found. Run 'tauso setup-genome --genome {genome}'")
    return gffutils.FeatureDB(paths["gff_db"])


def load_gff_pyranges(genome="GRCh38"):
    paths = get_paths(genome)
    if not os.path.exists(paths["gff"]):
        raise FileNotFoundError(f"GFF for {genome} not found. Run 'tauso setup-genome --genome {genome}'")
    return pr.read_gff3(paths["gff"])


def load_gtf_pyranges(genome="GRCh38"):
    paths = get_paths(genome)
    if not os.path.exists(paths["gtf"]):
        raise FileNotFoundError(f"GTF for {genome} not found. Run 'tauso setup-genome --genome {genome}'")
    return pr.read_gtf(paths["gtf"])


def _open_any_gtf(filepath: str, mode: str = "rt"):
    """Helper to transparently open either gzipped or uncompressed files."""
    if filepath.endswith(".gz"):
        return gzip.open(filepath, mode)
    return open(filepath, mode)


def _create_filtered_gtf(input_gtf_path: str, output_temp_path: str):
    """Streams a GTF and writes only header and 'gene' lines to the output path."""

    with _open_any_gtf(input_gtf_path, "rt") as f_in, open(output_temp_path, "wt") as f_out:
        for line in f_in:
            # Keep header lines
            if line.startswith("#"):
                f_out.write(line)
                continue

            # Fast split to check the Feature column (column 3)
            parts = line.split("\t", 3)
            if len(parts) > 2 and parts[2] == "gene":
                f_out.write(line)


def load_gtf_pyranges_gene_only(gtf_path: str):
    logger.info(f"[Load_GTF_assign_gene] Pre-filtering GTF only for genes {gtf_path}")

    # Generate a safe temporary file path
    fd, temp_path = tempfile.mkstemp(suffix=".gtf")
    os.close(fd)

    try:
        # Step 1: Filter the heavy GTF into the tiny temp file
        _create_filtered_gtf(gtf_path, temp_path)

        # Step 2: Load the tiny file into PyRanges
        logger.info(f"[Load_GTF_assign_gene] Loading filtered temporary GTF into RAM")
        gr = pr.read_gtf(temp_path)

    finally:
        # Step 3: Always clean up the temp file, even if the code crashes
        if os.path.exists(temp_path):
            os.remove(temp_path)

    logger.info(f"[Load_GTF_assign_gene] Filtered gene only regions, left with {len(gr)} rows in gr")
    return gr


def load_genome(genome="GRCh38"):
    """Returns the pyfaidx Fasta object."""
    paths = get_paths(genome)
    if not os.path.exists(paths["fasta"]):
        raise FileNotFoundError(f"FASTA for {genome} not found. Run 'tauso setup-genome --genome {genome}'")
    return Fasta(paths["fasta"])
