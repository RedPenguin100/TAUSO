"""Cytoplasmic rRNA reference targets for off-target scoring.

rRNA dominates total RNA, so RiboGreen / total-RNA-normalized assays are sensitive to
ASOs that hybridize rRNA — a signal the transcriptome off-target features miss (18S/28S
are absent from GRCh38; rRNA is ~zero-weighted by rRNA-depleted RNA-seq). The mature
cytoplasmic rRNA RefSeq sequences are fetched into the data dir by `tauso setup-rrna`.
"""

import logging
import urllib.request
from pathlib import Path

from tauso.data.data import get_data_dir
from tauso.genome.LocusInfo import LocusInfo

logger = logging.getLogger(__name__)

REFERENCE_FILENAME = "rrna_reference.fa"

# feature name -> (RefSeq accession, expected length for a sanity check)
RRNA_ACCESSIONS = {
    "rRNA_18S": ("NR_003286.4", 1869),
    "rRNA_5_8S": ("NR_003285.3", 157),
    "rRNA_28S": ("NR_003287.4", 5070),
    "rRNA_5S": ("NR_023363.1", 119),
}

_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype=fasta&retmode=text"


def reference_path() -> Path:
    return Path(get_data_dir()) / REFERENCE_FILENAME


def fetch_rrna_reference(path: Path | None = None, overwrite: bool = False) -> Path:
    """Download the cytoplasmic rRNA RefSeq sequences into a single FASTA (needs network)."""
    path = path or reference_path()
    if path.exists() and not overwrite:
        return path

    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    for name, (acc, expected_len) in RRNA_ACCESSIONS.items():
        logger.info("Fetching %s (%s) from NCBI...", name, acc)
        with urllib.request.urlopen(_EFETCH_URL.format(acc=acc), timeout=60) as resp:
            raw = resp.read().decode("ascii")
        seq = "".join(l.strip() for l in raw.splitlines() if l and not l.startswith(">"))
        if not seq:
            raise RuntimeError(f"Empty sequence fetched for {name} ({acc}).")
        if abs(len(seq) - expected_len) > 5:
            logger.warning("%s (%s) length %d differs from expected ~%d nt.", name, acc, len(seq), expected_len)
        lines.append(f">{name} {acc}")
        lines.extend(seq[i : i + 70] for i in range(0, len(seq), 70))
    path.write_text("\n".join(lines) + "\n")
    logger.info("Wrote %d rRNA sequences to %s.", len(RRNA_ACCESSIONS), path)
    return path


def load_rrna_sequences(path: Path | None = None) -> dict[str, str]:
    """Parse the rRNA reference FASTA into {feature_name: sequence}."""
    path = path or reference_path()
    if not path.exists():
        raise FileNotFoundError(f"rRNA reference FASTA not found at {path}. Run `tauso setup-rrna`.")
    seqs: dict[str, str] = {}
    name = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(chunks)
            name = line[1:].split()[0]
            chunks = []
        elif line.strip():
            chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def get_rrna_loci(path: Path | None = None) -> dict[str, LocusInfo]:
    """Return {feature_name: LocusInfo} for each rRNA species (LocusInfo.full_mrna = sequence)."""
    return {name: LocusInfo(seq=seq) for name, seq in load_rrna_sequences(path).items()}
