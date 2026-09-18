"""How the pipeline stores text in a DataFrame, and what it hands back when it is done."""

import numpy as np
import pandas as pd
import pyarrow as pa

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


def release_arrow_memory():
    """Give pyarrow's decoded pages back to the operating system.

    Reading a parquet file leaves what it decoded in pyarrow's own pool, ready for the next
    read. Those pages are not Python objects, so neither the garbage collector nor glibc's
    allocator can return them: after the half-life table is read, 371 MB sits there for the
    rest of the run. Whoever has finished reading has no use for them.
    """
    pa.default_memory_pool().release_unused()
