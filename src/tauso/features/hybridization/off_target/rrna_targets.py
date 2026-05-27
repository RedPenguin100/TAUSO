"""Cytoplasmic rRNA reference targets for off-target scoring.

rRNA dominates total RNA, so RiboGreen / total-RNA-normalized assays are sensitive to
ASOs that hybridize rRNA — a signal the transcriptome off-target features miss (18S/28S
are absent from GRCh38; rRNA is ~zero-weighted by rRNA-depleted RNA-seq). The mature
cytoplasmic rRNA RefSeq sequences are fetched into the data dir by `tauso setup-rrna`.

FASTA download/parse/write live in tauso.data.ncbi and tauso.genome.fasta; this module
only holds the rRNA-specific accessions and the loci injected into off-target scoring.
"""

import logging
from io import StringIO
from pathlib import Path

from Bio import SeqIO

from ....data.data import get_data_dir
from ....data.ncbi import fetch_nuccore_fasta
from ....genome.fasta import read_fasta, write_fasta
from ....genome.LocusInfo import LocusInfo

logger = logging.getLogger(__name__)

REFERENCE_FILENAME = "rrna_reference.fa"

# feature name -> (RefSeq accession, expected length for a sanity check)
RRNA_ACCESSIONS = {
    "rRNA_18S": ("NR_003286.4", 1869),
    "rRNA_5_8S": ("NR_003285.3", 157),
    "rRNA_28S": ("NR_003287.4", 5070),
    "rRNA_5S": ("NR_023363.1", 119),
}


def reference_path() -> Path:
    return Path(get_data_dir()) / REFERENCE_FILENAME


def fetch_rrna_reference(path: Path | None = None, overwrite: bool = False) -> Path:
    """Download the cytoplasmic rRNA RefSeq sequences into a single FASTA (needs network)."""
    path = path or reference_path()
    if path.exists() and not overwrite:
        return path

    sequences: dict[str, str] = {}
    for name, (acc, expected_len) in RRNA_ACCESSIONS.items():
        logger.info("Fetching %s (%s) from NCBI...", name, acc)
        seq = str(SeqIO.read(StringIO(fetch_nuccore_fasta(acc)), "fasta").seq)
        if not seq:
            raise RuntimeError(f"Empty sequence fetched for {name} ({acc}).")
        if abs(len(seq) - expected_len) > 5:
            logger.warning("%s (%s) length %d differs from expected ~%d nt.", name, acc, len(seq), expected_len)
        sequences[name] = seq

    write_fasta(sequences, path)
    logger.info("Wrote %d rRNA sequences to %s.", len(sequences), path)
    return path


def load_rrna_sequences(path: Path | None = None) -> dict[str, str]:
    """Read the rRNA reference FASTA into {feature_name: sequence}."""
    path = path or reference_path()
    if not path.exists():
        raise FileNotFoundError(f"rRNA reference FASTA not found at {path}. Run `tauso setup-rrna`.")
    return read_fasta(path)


def get_rrna_loci(path: Path | None = None) -> dict[str, LocusInfo]:
    """Return {feature_name: LocusInfo} for each rRNA species (LocusInfo.full_mrna = sequence)."""
    return {name: LocusInfo(seq=seq) for name, seq in load_rrna_sequences(path).items()}
