"""Every hit of one search, for tests that compare runs against each other.

Production reduces hits as they arrive and never holds them all; these tests
want them all, so this is here rather than beside the code under test.
"""

import pandas as pd
import pyrisearch_tauso

from tauso.features.hybridization.fast_hybridization import (
    ASO_TARGET_MATRIX,
    EXTENSION_PENALTY,
    RISEARCH_COLUMNS,
)


def risearch_hits_dataframe(
    query_id_seq_pairs,
    target_file_path,
    *,
    matrix=ASO_TARGET_MATRIX,
    minimum_score=900,
    neighborhood=0,
    transpose=False,
):
    if not query_id_seq_pairs:
        return pd.DataFrame(columns=list(RISEARCH_COLUMNS))

    return pyrisearch_tauso.hits_table(
        queries=query_id_seq_pairs,
        targets=target_file_path,
        min_score=minimum_score,
        matrix=matrix,
        extension_penalty=EXTENSION_PENALTY,
        neighborhood=neighborhood,
        transpose=transpose,
    ).to_pandas()
