import pandas as pd

from notebooks.data.OligoAI.parse_chemistry import assign_chemistry
from notebooks.data.utility import standardize_cell_line
from notebooks.preprocessing import process_oligo_data_rename
from tauso.data.consts import SEQUENCE, CELL_LINE, CANONICAL_GENE, CELL_LINE_ORGANISM, INHIBITION, VOLUME

def standardize_oligo_ai_data(df: pd.DataFrame) -> pd.DataFrame:

    df = process_oligo_data_rename(df)

    rename_scheme = {'aso_sequence_5_to_3': SEQUENCE, 'Canonical Gene Name': CANONICAL_GENE, 'cell_line': CELL_LINE,
                     'cell_line_species': CELL_LINE_ORGANISM, 'inhibition_percent': INHIBITION, 'dosage': VOLUME}
    df = df.rename(columns=rename_scheme)

    df = df[df['steric_blocking'] == False]
    df = df[df[INHIBITION].notna()]
    df = standardize_cell_line(df)
    df = assign_chemistry(df)

    return df