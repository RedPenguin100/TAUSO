import pandas as pd
import pytest

from notebooks.features.calculate_features import get_partition_dir, merge_partitions, partition_data
from tauso.data.consts import CANONICAL_GENE, PS_PATTERN
from tauso.populate.calculators.calculator import Calculator

pytestmark = pytest.mark.integration


@pytest.fixture
def small_df():
    """Two genes, 5 rows each — enough to test partition + merge without any external data."""
    rows = []
    for i, gene in enumerate(["KRAS", "BRAF"]):
        for j in range(5):
            rows.append({
                "index_oligo": i * 5 + j + 1,
                CANONICAL_GENE: gene,
                PS_PATTERN: "***ddd***",
            })
    return pd.DataFrame(rows)


def test_partition_data_is_disjoint_and_complete(small_df):
    p0 = partition_data(small_df, k=0, n=2)
    p1 = partition_data(small_df, k=1, n=2)

    genes_p0 = set(p0[CANONICAL_GENE].unique())
    genes_p1 = set(p1[CANONICAL_GENE].unique())

    assert genes_p0 & genes_p1 == set(), "Partitions must not share genes"
    assert genes_p0 | genes_p1 == {"KRAS", "BRAF"}, "Partitions must cover all genes"
    assert len(p0) + len(p1) == len(small_df), "No rows lost or duplicated"


def test_partition_and_merge_backbone_features(small_df, tmp_path):
    version = "oligo"
    main_dir = tmp_path / f"saved_features_{version}"
    index_col = f"index_{version}"

    for k in range(2):
        part_data = partition_data(small_df, k=k, n=2)
        part_dir = get_partition_dir(main_dir, k=k, n=2)

        calc = Calculator(
            data=part_data,
            data_version=version,
            overwrite=False,
            cpus=1,
            get_feature_dir=lambda v, d=part_dir: d,
        )
        calc.calculate_backbone_features()

        assert (part_dir / "ps_po_percentage.csv").exists(), f"Partition {k} should have saved ps_po_percentage.csv"
        assert len(pd.read_csv(part_dir / "ps_po_percentage.csv")) == 5

    merge_partitions(main_dir, n_partitions=2, index_col=index_col)

    for feature in ["ps_po_percentage", "ps_end_score", "ps_max_consecutive_po"]:
        merged = pd.read_csv(main_dir / f"{feature}.csv")
        assert len(merged) == len(small_df), f"{feature}: expected {len(small_df)} rows after merge"
        assert set(merged[index_col]) == set(small_df[index_col]), f"{feature}: merged indices must match original"
