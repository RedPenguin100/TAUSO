import pandas as pd

from notebooks.data.OligoAI.parse_chemistry import assign_chemistry
from notebooks.data.utility import standardize_cell_line
from notebooks.preprocessing import process_oligo_data_rename
from tauso.data.consts import CANONICAL_GENE, INHIBITION

def standardize_oligo_ai_data(df: pd.DataFrame) -> pd.DataFrame:
    df = process_oligo_data_rename(df)
    df = df.rename(columns={'Canonical Gene Name': CANONICAL_GENE})
    df = df[df['steric_blocking'] == False]
    df = df[df[INHIBITION].notna()]
    df = standardize_cell_line(df)
    df = assign_chemistry(df)
    return df
