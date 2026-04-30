from pathlib import Path


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
    rna_sequence = "".join(sequence_parts).upper().replace("T", "U")

    if not rna_sequence:
        raise ValueError(f"Validation Failed: Sequence for gene '{gene_name}' is empty.")

    return gene_name, rna_sequence
