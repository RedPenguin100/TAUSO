import pandas as pd

from tauso.data.consts import (
    CELL_LINE_DEPMAP,
    CELL_LINE_DEPMAP_PROXY,
    CELL_LINE,
    CELL_LINE_TO_DEPMAP_PROXY_DICT,
    CELL_LINE_TO_DEPMAP,
    CANONICAL_GENE,
    SEQUENCE,
    SENSE_START,
)


def standardize_cell_line(df: pd.DataFrame) -> pd.DataFrame:
    result = df.copy()
    result[CELL_LINE_DEPMAP_PROXY] = result[CELL_LINE].map(CELL_LINE_TO_DEPMAP_PROXY_DICT)
    result[CELL_LINE_DEPMAP] = result[CELL_LINE_DEPMAP_PROXY].map(CELL_LINE_TO_DEPMAP)

    return result


import pandas as pd

# Translation table: converts A->U, C->G, G->C, T/U->A
COMP_U_TABLE = bytes.maketrans(b"ACGTUacgtu", b"UGCAaugcaa")

def get_antisense_u(seq: str) -> str:
    """Returns the reverse complement RNA (sense) sequence."""
    return seq.encode().translate(COMP_U_TABLE)[::-1].decode()

def find_sense_start(row, gene_to_data):
    """Finds the 0-indexed start position of the sense sequence in the mRNA."""
    gene = row[CANONICAL_GENE]
    antisense_seq = row[SEQUENCE]

    # Return -1 if gene isn't in our dictionary
    if pd.isna(gene) or gene not in gene_to_data:
        return -1

    # Get the mRNA and ensure it's formatted as uppercase RNA
    mrna = gene_to_data[gene].full_mrna.upper().replace('T', 'U')

    # Convert the antisense input into its sense counterpart
    sense_seq = get_antisense_u(antisense_seq)

    # .find() natively returns the starting index, or -1 if not found
    return mrna.find(sense_seq)

def get_populated_df_with_sense_start(df, gene_to_data):
    """Main wrapper to apply the logic to the DataFrame."""
    # Work on a copy to prevent SettingWithCopy warnings
    result_df = df.copy()

    result_df[SENSE_START] = result_df.apply(
        lambda row: find_sense_start(row, gene_to_data),
        axis=1
    )

    return result_df