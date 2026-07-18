import os

import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Implicit n_jobs helper
# Tests call get_n_jobs() inline — no fixture parameter needed.
# Value is set once at session start from the --n-jobs CLI option.
# ---------------------------------------------------------------------------
_n_jobs: list[int] = [1]


def get_n_jobs() -> int:
    return _n_jobs[0]


@pytest.fixture(scope="session", autouse=True)
def _configure_n_jobs(request):
    _n_jobs[0] = request.config.getoption("--n-jobs")


from notebooks.consts import OLIGO_CSV_INDEXED
from notebooks.data.OligoAI.parse_chemistry import assign_chemistry
from notebooks.data.OligoAI.utility import standardize_oligo_ai_data

from tauso.algorithms.genomic_context_windows import (
    add_external_mrna_and_context_columns,
)
from tauso.data.consts import *
from tauso.data.consts import CELL_LINE_DEPMAP
from tauso.data.data import get_data_dir
from tauso.features.codon_usage.find_cai_reference import load_cell_line_gene_expression
from tauso.features.context.mrna_halflife import HalfLifeProvider, load_halflife_mapping
from tauso.features.rbp.load_rbp import load_attract_data
from tauso.genome.read_human_genome import get_locus_to_data_dict
from tauso.genome.TranscriptMapper import build_gene_sequence_registry
from tauso.populate.populate_structure import get_populated_df_with_structure_features
from tauso.timer import Timer


@pytest.fixture(scope="session")
def transcriptomes(base_data, target_genes):
    """Loads cell line gene expression data."""
    with Timer("Load Transcriptomes"):
        data_dir = get_data_dir()
        expression_dir = os.path.join(data_dir, "processed_expression")
        depmap_ids = list(set(base_data[CELL_LINE_DEPMAP]))

        trans = load_cell_line_gene_expression(depmap_ids, target_genes, expression_dir=expression_dir)
        return trans


@pytest.fixture(scope="session")
def rbp_assets():
    """Loads ATTRACT RBP data returning (rbp_map, pwm_db)."""
    with Timer("Load ATTRACT RBP data"):
        return load_attract_data()


FLANK_SIZES_PREMRNA = [20, 30, 40, 50, 60, 70]
CDS_WINDOWS = [20, 30, 40, 50, 60, 70]


@pytest.fixture(scope="session")
def raw_oligo_data():
    with Timer("Load CSV Data"):
        data = pd.read_csv(OLIGO_CSV_INDEXED, engine="pyarrow")

    if "index_oligo" not in data.columns:
        raise ValueError("index_oligo column not found")
    return data


@pytest.fixture(scope="session")
def base_data(raw_oligo_data):
    """Handles basic renaming, filtering, cell line mapping, and one-hot encoding."""
    with Timer("Rename Columns and Basic Filtering"):
        data = standardize_oligo_ai_data(raw_oligo_data)

    with Timer("Transfection One Hot Encoding"):
        from tauso.populate.populate_context import populate_transfection

        data, _ = populate_transfection(data)

    return data


@pytest.fixture(scope="session")
def target_genes(base_data):
    """Extracts unique target genes from the base data."""
    return base_data[CANONICAL_GENE_NAME].unique().tolist()


@pytest.fixture(scope="session")
def gene_to_data(target_genes):
    """Fetches the locus to data dictionary based on target genes."""
    with Timer("Get Locus to Data Dict"):
        return get_locus_to_data_dict(include_introns=True, gene_subset=target_genes)


@pytest.fixture(scope="session")
def gene_to_data_full(target_genes):
    """Fetches the locus to data dictionary for all genes"""
    with Timer("Get Locus to Data Dict"):
        return get_locus_to_data_dict(include_introns=True)


@pytest.fixture(scope="session")
def half_life_provider():
    with Timer("Get mRNA Half Life provider"):
        mapping = load_halflife_mapping()
        provider = HalfLifeProvider(mapping)
        return provider


@pytest.fixture(scope="session")
def structure_data(base_data, target_genes, gene_to_data):
    """Populates the dataframe with structure features."""
    with Timer("Populate DF with Structure Features"):
        return get_populated_df_with_structure_features(base_data, target_genes, gene_to_data)


@pytest.fixture(scope="session")
def chemistry_data(structure_data):
    return assign_chemistry(structure_data)


@pytest.fixture(scope="session")
def ref_registry(target_genes, gene_to_data):
    """Builds the Gene Sequence Registry."""
    with Timer("Build Gene Sequence Registry"):
        return build_gene_sequence_registry(genes=target_genes, gene_to_data=gene_to_data)


@pytest.fixture(scope="session")
def final_data(structure_data, ref_registry):
    """Runs the optimized context generator to add external mRNA and context columns."""
    with Timer("Add External mRNA & Context Columns"):
        return add_external_mrna_and_context_columns(
            df=structure_data,
            gene_registry=ref_registry,
            flank_sizes_premrna=FLANK_SIZES_PREMRNA,
            flank_sizes_cds=CDS_WINDOWS,
        )


@pytest.fixture(scope="session")
def mini_sampled_data(request, final_data):
    """Samples final_data once per (n_samples, session) and caches the result.

    Tests must call .copy() before mutating the returned DataFrame.
    """
    n_samples = getattr(request, "param", 1000)
    actual_samples = min(n_samples, len(final_data))
    return final_data.sample(n=actual_samples, random_state=42).copy()


@pytest.fixture(scope="session")
def common_gene_sample(request, base_data):
    """A small, deterministic slice: `max_per_gene` ASOs from each of the `n_genes` most-
    represented genes. Parametrize indirectly via (n_genes, max_per_gene); defaults to (2, 6).

    Tests must call .copy() before mutating the returned DataFrame.
    """
    n_genes, max_per_gene = getattr(request, "param", (2, 6))
    genes = base_data[CANONICAL_GENE_NAME].value_counts().index[:n_genes].tolist()
    sub = base_data[base_data[CANONICAL_GENE_NAME].isin(genes)]
    return sub.groupby(CANONICAL_GENE_NAME, group_keys=False).head(max_per_gene).reset_index(drop=True)


@pytest.fixture
def sampled_base_data(request, base_data):
    """Samples the base data right before the test runs."""
    n_samples = getattr(request, "param", 1000)
    actual_samples = min(n_samples, len(base_data))
    return base_data.sample(n=actual_samples, random_state=42).copy()
