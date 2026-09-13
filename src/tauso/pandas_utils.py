"""Frame helpers shared by the feature code."""

import pandas as pd


def add_columns(df: pd.DataFrame, columns: dict) -> pd.DataFrame:
    """Return `df` with `columns` ({name: values}) attached.

    A name already on the frame has its values replaced and keeps its position; the rest are
    appended in the order given. The frame is written in place, so it never has to be copied
    however wide it has grown.
    """
    for name, values in columns.items():
        df[name] = values
    return df
