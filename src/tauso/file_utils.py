import gzip
import warnings
from pathlib import Path

from Bio import SeqIO

from .data.data import get_paths
from .timer import Timer


def get_fasta_dict_from_path(fasta_path: Path):
    with Timer() as timer:
        if fasta_path.suffix == ".gz":
            warnings.warn(
                f"Fasta is compressed, consider decompressing for performance. To unzip, run gunzip {fasta_path}",
                stacklevel=2,
            )
            with gzip.open(str(fasta_path), "rt") as handle:
                fasta_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        else:
            with open(str(fasta_path), "r") as handle:
                fasta_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    print(f"Time took to read fasta: {timer.elapsed_time}")
    return fasta_dict


def read_genome_fasta_dict(genome: str):
    paths = get_paths(genome)
    genome_fasta = Path(paths["fasta"])
    if genome_fasta.is_file():
        return get_fasta_dict_from_path(genome_fasta)

    raise FileNotFoundError(f"Did not find {genome_fasta}, please consider the README.md")


def read_human_genome_fasta_dict():
    return read_genome_fasta_dict(genome="GRCh38")
