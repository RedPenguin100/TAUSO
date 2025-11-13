import warnings
import gzip

from pathlib import Path
from Bio import SeqIO

from .consts import HUMAN_TRANSCRIPTS_FASTA, HUMAN_TRANSCRIPTS_FASTA_GZ, HUMAN_GENOME_FASTA, \
    HUMAN_GENOME_FASTA_GZ, YEAST_FASTA_PATH
from .timer import Timer


def get_fasta_dict_from_path(fasta_path: Path):
    with Timer() as timer:
        if fasta_path.suffix == ".gz":
            warnings.warn(
                f"Fasta is compressed, consider decompressing for performance. To unzip, run gunzip {fasta_path}")
            with gzip.open(str(fasta_path), 'rt') as handle:
                fasta_dict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
        else:
            with open(str(fasta_path), 'r') as handle:
                fasta_dict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
    print(f"Time took to read fasta: {timer.elapsed_time}")
    return fasta_dict


def read_human_transcriptome_fasta_dict():
    if HUMAN_TRANSCRIPTS_FASTA.is_file():
        return get_fasta_dict_from_path(HUMAN_TRANSCRIPTS_FASTA)

    if HUMAN_TRANSCRIPTS_FASTA_GZ.is_file():
        return get_fasta_dict_from_path(HUMAN_TRANSCRIPTS_FASTA_GZ)

    raise FileNotFoundError(
        f"Did not find {HUMAN_TRANSCRIPTS_FASTA} or {HUMAN_TRANSCRIPTS_FASTA_GZ}, please consider the README.md")


def read_human_genome_fasta_dict():
    if HUMAN_GENOME_FASTA.is_file():
        return get_fasta_dict_from_path(HUMAN_GENOME_FASTA)

    if HUMAN_GENOME_FASTA_GZ.is_file():
        return get_fasta_dict_from_path(HUMAN_GENOME_FASTA_GZ)

    raise FileNotFoundError(
        f"Did not find {HUMAN_GENOME_FASTA} or {HUMAN_GENOME_FASTA_GZ}, please consider the README.md")


def read_yeast_genome_fasta_dict():
    if YEAST_FASTA_PATH.is_file():
        return get_fasta_dict_from_path(YEAST_FASTA_PATH)
    raise FileNotFoundError(f"Did not find {YEAST_FASTA_PATH}, please consider the README.md")
