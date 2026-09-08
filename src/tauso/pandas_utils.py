"""Frame helpers shared by the feature code."""

import pandas as pd


def add_columns(df: pd.DataFrame, columns: dict) -> pd.DataFrame:
    """Return `df` with `columns` ({name: values}) attached in one pass.

    A name already on the frame has its values replaced and keeps its position; the rest are
    appended together. Assigning a column at a time inserts it into the existing frame, and
    pandas rebuilds the frame on every insert.
    """
    fresh = {name: values for name, values in columns.items() if name not in df.columns}
    for name, values in columns.items():
        if name not in fresh:
            df[name] = values
    if fresh:
        df = pd.concat([df, pd.DataFrame(fresh, index=df.index)], axis=1)
    return df
