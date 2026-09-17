"""How the pipeline stores text in a DataFrame."""

import numpy as np
import pandas as pd

# A column of Python strings costs about 60 bytes of header per cell on top of the text;
# Arrow keeps one buffer per column. This is the variant that leaves missing values as NaN
# and returns numpy arrays from the string methods, so a column reads as it always did.
ARROW_STRINGS = pd.StringDtype("pyarrow", na_value=np.nan)


def arrow_strings(df):
    """Store every text column of `df` as Arrow strings rather than Python objects.

    Only columns holding nothing but strings (and missing values) are converted; the rest
    stay as they are.
    """
    for column in df.columns[(df.dtypes == object).to_numpy()]:
        if pd.api.types.infer_dtype(df[column], skipna=True) == "string":
            df[column] = df[column].astype(ARROW_STRINGS)
    return df
