import logging
import time

import pytest

from tauso.features.context.ribo_seq import add_genomic_coordinates
from tauso.populate.populate_context import populate_ribo_seq


@pytest.mark.parametrize("mini_sampled_data", [500], indirect=True)
def test_ribo_seq(mini_sampled_data, mapper, dataframe_regression):
    """Regression test: ribo-seq feature values must match the stored baseline."""
    data = mini_sampled_data.copy()

    t0 = time.perf_counter()
    data = add_genomic_coordinates(data, mapper)
    t_coords = time.perf_counter() - t0

    t1 = time.perf_counter()
    data, feature_cols = populate_ribo_seq("human", data)
    t_ribo = time.perf_counter() - t1

    logging.getLogger(__name__).info(
        "add_genomic_coordinates: %.2fs | populate_ribo_seq: %.2fs | total: %.2fs",
        t_coords,
        t_ribo,
        t_coords + t_ribo,
    )

    check_cols = ["index_oligo"] + feature_cols
    dataframe_regression.check(data[check_cols])
