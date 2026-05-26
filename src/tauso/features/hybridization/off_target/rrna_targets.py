"""Cytoplasmic rRNA reference targets for off-target ("denominator distortion") scoring.

RiboGreen-based ASO potency assays normalize the target mRNA signal to *total RNA*,
which is ~85% rRNA by mass. An ASO with RNase-H-competent complementarity to rRNA
therefore distorts the assay denominator (and/or sequesters RNase H1) — an
experimental/competition artifact rather than a true property of target engagement.

The standard transcriptome off-target features cannot see this:
  * 18S and 28S (the bulk of rRNA mass) are ABSENT from GRCh38 — the 45S rDNA repeat
    sits in the unassembled acrocentric NORs (chr13/14/15/21/22 short arms), so they
    are not in the gene annotation at all.
  * Even the 5S / 5.8S that ARE annotated get ~zero weight, because the off-target
    features weight binding by DepMap RNA-seq TPM and RNA-seq libraries are
    rRNA-depleted (poly-A selected / Ribo-Zero), so rRNA reads ~0 TPM.

So we supply the mature cytoplasmic rRNA sequences explicitly from RefSeq and score
each ASO against them with the existing single-gene RIsearch path
(``off_target_specific_seq_pandarallel``), producing features
``off_target_single_{name}_c{cutoff}``.

The reference FASTA is fetched once into the data dir. Regenerate with::

    python -m tauso.features.hybridization.off_target.rrna_targets
"""

import logging
import urllib.request
from pathlib import Path

from tauso.data.data import get_data_dir
from tauso.genome.LocusInfo import LocusInfo

logger = logging.getLogger(__name__)

REFERENCE_FILENAME = "rrna_reference.fa"

# Mature cytoplasmic rRNA species, keyed by the name used for the feature
# (off_target_single_<name>_c<cutoff>) -> curated human RefSeq accession.
# Lengths are recorded for a post-fetch sanity check.
RRNA_ACCESSIONS = {
    "rRNA_18S": ("NR_003286.4", 1869),   # RNA18SN5
    "rRNA_5_8S": ("NR_003285.3", 157),   # RNA5-8SN5
    "rRNA_28S": ("NR_003287.4", 5070),   # RNA28SN5
    "rRNA_5S": ("NR_023363.1", 119),     # RNA5S1
}

_EFETCH_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=nuccore&id={acc}&rettype=fasta&retmode=text"
)


def reference_path() -> Path:
    """Path to the rRNA reference FASTA inside the data dir."""
    return Path(get_data_dir()) / REFERENCE_FILENAME


def fetch_rrna_reference(path: Path | None = None, overwrite: bool = False) -> Path:
    """Download the cytoplasmic rRNA RefSeq sequences into a single FASTA.

    Each record is written with a ``>{name} {accession}`` header so the loader can
    key on the feature name. Network access (NCBI E-utilities) is required.
    """
    path = path or reference_path()
    if path.exists() and not overwrite:
        logger.info("rRNA reference already present at %s (use overwrite=True to refetch).", path)
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
        raise FileNotFoundError(
            f"rRNA reference FASTA not found at {path}. "
            "Run `python -m tauso.features.hybridization.off_target.rrna_targets` to fetch it."
        )
    seqs: dict[str, str] = {}
    name = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(chunks)
            name = line[1:].split()[0]  # first token = feature name
            chunks = []
        elif line.strip():
            chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def get_rrna_loci(path: Path | None = None) -> dict[str, LocusInfo]:
    """Return {feature_name: LocusInfo} for each cytoplasmic rRNA species.

    ``LocusInfo(seq=...)`` sets ``full_mrna`` directly, which is what the RIsearch
    off-target scorer (``_apply_risearch_scoring``) reads off each target.
    """
    return {name: LocusInfo(seq=seq) for name, seq in load_rrna_sequences(path).items()}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    fetch_rrna_reference(overwrite=True)
