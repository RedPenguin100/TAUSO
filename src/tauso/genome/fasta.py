from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..util import dna_to_rna


def read_fasta(file_path: str | Path) -> dict[str, str]:
    """Parse a (multi-record) FASTA into {record_id: sequence}.

    record_id is the first whitespace-delimited token of each header line.
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found at: {path}")
    return {record.id: str(record.seq) for record in SeqIO.parse(str(path), "fasta")}


def write_fasta(records: dict[str, str], file_path: str | Path) -> None:
    """Write {record_id: sequence} to a FASTA file (creating parent dirs)."""
    path = Path(file_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    seq_records = [SeqRecord(Seq(seq), id=name, description="") for name, seq in records.items()]
    SeqIO.write(seq_records, str(path), "fasta")


def read_single_rna_fasta(file_path: str | Path) -> tuple[str, str]:
    """
    Reads a FASTA file, verifies it contains exactly one sequence.
    Returns: (gene_name, rna_sequence)
    Raises: ValueError if format is invalid or multiple genes exist.
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found at: {path}")

    gene_name = ""
    sequence_parts = []
    header_count = 0

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue  # Skip empty lines

            if line.startswith(">"):
                header_count += 1
                if header_count > 1:
                    raise ValueError("Validation Failed: Multiple sequences found in file.")

                # Extract name: remove '>' and split by space to drop extra descriptions
                full_header = line[1:].strip()
                if not full_header:
                    raise ValueError("Validation Failed: Gene header exists but has no name.")

                gene_name = full_header.split()[0]
            else:
                # Catch malformed files that have sequence text before the > header
                if header_count == 0:
                    raise ValueError("Validation Failed: Sequence data found before a header line.")
                sequence_parts.append(line)

    if header_count == 0:
        raise ValueError("Validation Failed: No FASTA header (starting with '>') found.")

    # Process sequence: Uppercase and T -> U
    rna_sequence = dna_to_rna("".join(sequence_parts))

    if not rna_sequence:
        raise ValueError(f"Validation Failed: Sequence for gene '{gene_name}' is empty.")

    return gene_name, rna_sequence
