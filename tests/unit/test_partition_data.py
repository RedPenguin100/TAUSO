import pandas as pd
from notebooks.features.calculate_features import partition_data

from tauso.data.consts import CANONICAL_GENE_NAME


def test_partition_data_is_disjoint_and_complete():
    """A gene must land in exactly one partition, or a cluster run loses or duplicates rows."""
    df = pd.DataFrame([{CANONICAL_GENE_NAME: gene} for gene in ["KRAS", "BRAF"] for _ in range(5)])

    p0 = partition_data(df, k=0, n=2)
    p1 = partition_data(df, k=1, n=2)

    genes_p0 = set(p0[CANONICAL_GENE_NAME].unique())
    genes_p1 = set(p1[CANONICAL_GENE_NAME].unique())

    assert genes_p0 & genes_p1 == set(), "Partitions must not share genes"
    assert genes_p0 | genes_p1 == {"KRAS", "BRAF"}, "Partitions must cover all genes"
    assert len(p0) + len(p1) == len(df), "No rows lost or duplicated"
