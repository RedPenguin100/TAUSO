import os

import pandas as pd
import pytest

from notebooks.consts import OLIGO_CSV_INDEXED
from tauso.algorithms.genomic_context_windows import (
    add_external_mrna_and_context_columns,
)
from tauso.data.consts import *
from tauso.data.consts import CELL_LINE_DEPMAP
from tauso.data.data import get_data_dir, get_paths
from tauso.features.codon_usage.find_cai_reference import load_cell_line_gene_expression
from tauso.features.context.mrna_halflife import HalfLifeProvider, load_halflife_mapping
from tauso.features.rbp.load_rbp import load_attract_data
from tauso.features.rbp.pwm_helper import build_rbp_expression_matrix
from tauso.genome.read_human_genome import get_locus_to_data_dict
from tauso.genome.TranscriptMapper import (
    GeneCoordinateMapper,
    build_gene_sequence_registry,
)
from tauso.populate.populate_structure import get_populated_df_with_structure_features
from tauso.timer import Timer


@pytest.fixture(scope="session")
def transcriptomes(base_data, target_genes):
    """Loads cell line gene expression data."""
    with Timer("Load Transcriptomes"):
        data_dir = get_data_dir()
        expression_dir = os.path.join(data_dir, "processed_expression")
        depmap_ids = list(set(base_data[CELL_LINE_DEPMAP]))

        trans = load_cell_line_gene_expression(
            depmap_ids, target_genes, expression_dir=expression_dir
        )
        return trans


@pytest.fixture(scope="session")
def rbp_assets():
    """Loads ATTRACT RBP data returning (rbp_map, pwm_db)."""
    with Timer("Load ATTRACT RBP data"):
        return load_attract_data()


@pytest.fixture(scope="session")
def expression_matrix(base_data, rbp_assets):
    """Builds the RBP expression matrix."""
    with Timer("Build RBP Expression Matrix"):
        _, pwm_db = rbp_assets
        return build_rbp_expression_matrix(df=base_data, pwm_db=pwm_db)


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
    with Timer("Rename Columns & Basic Filtering"):
        rename_scheme = {
            "aso_sequence_5_to_3": SEQUENCE,
            "Canonical Gene Name": CANONICAL_GENE,
            "cell_line": CELL_LINE,
            "cell_line_species": CELL_LINE_ORGANISM,
            "inhibition_percent": INHIBITION,
            "dosage": VOLUME,
        }
        data = raw_oligo_data.rename(columns=rename_scheme)

        data = data[data["steric_blocking"] == False]
        data = data[~data[CANONICAL_GENE].str.contains(";")]
        data = data[data[INHIBITION].notna()]

    with Timer("Cell Line Mapping & One-Hot Encoding"):
        data[CELL_LINE_DEPMAP_PROXY] = data[CELL_LINE].map(
            CELL_LINE_TO_DEPMAP_PROXY_DICT
        )
        data[CELL_LINE_DEPMAP] = data[CELL_LINE_DEPMAP_PROXY].map(CELL_LINE_TO_DEPMAP)

        binary_features = pd.get_dummies(data["transfection_method"])
        data = pd.concat([data, binary_features], axis=1)

    return data


@pytest.fixture(scope="session")
def target_genes(base_data):
    """Extracts unique target genes from the base data."""
    return base_data[CANONICAL_GENE].unique().tolist()


@pytest.fixture(scope="session")
def mapper():
    """Initializes and returns the GeneCoordinateMapper."""
    with Timer("Path & Mapper Initialization"):
        paths = get_paths("GRCh38")
        return GeneCoordinateMapper(paths["db"])


@pytest.fixture(scope="session")
def gene_to_data(target_genes):
    """Fetches the locus to data dictionary based on target genes."""
    with Timer("Get Locus to Data Dict"):
        return get_locus_to_data_dict(include_introns=True, gene_subset=target_genes)


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
        return get_populated_df_with_structure_features(
            base_data, target_genes, gene_to_data
        )


@pytest.fixture(scope="session")
def ref_registry(target_genes, gene_to_data, mapper):
    """Builds the Gene Sequence Registry."""
    with Timer("Build Gene Sequence Registry"):
        return build_gene_sequence_registry(
            genes=target_genes, gene_to_data=gene_to_data, mapper=mapper
        )


@pytest.fixture(scope="session")
def final_data(structure_data, mapper, ref_registry):
    """Runs the optimized context generator to add external mRNA and context columns."""
    with Timer("Add External mRNA & Context Columns"):
        return add_external_mrna_and_context_columns(
            df=structure_data,
            mapper=mapper,
            gene_registry=ref_registry,
            flank_sizes_premrna=FLANK_SIZES_PREMRNA,
            flank_sizes_cds=CDS_WINDOWS,
        )


@pytest.fixture
def mini_sampled_data(request, final_data):
    """Samples the fully processed DataFrame right before the test runs."""
    n_samples = getattr(request, "param", 1000)
    actual_samples = min(n_samples, len(final_data))
    return final_data.sample(n=actual_samples, random_state=42).copy()
