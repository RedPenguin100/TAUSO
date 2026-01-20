import pytest
import pandas as pd

# Imports from your project
from notebooks.consts import UPDATED_CSV
from notebooks.preprocessing import preprocess_aso_data, get_unique_genes
from tauso.data.data import get_paths
from tauso.genome.TranscriptMapper import GeneCoordinateMapper, build_gene_sequence_registry
from tauso.genome.read_human_genome import get_locus_to_data_dict
from tauso.algorithms.genomic_context_windows import add_external_mrna_and_context_columns


def test_pipeline_data_regression(data_regression):
        # 1. Load & Preprocess
    full_data = preprocess_aso_data(UPDATED_CSV)
    df_subset = full_data.copy()

    # 2. PREPARATION (Outside the function)
    # A. Get unique genes involved
    genes_subset = get_unique_genes(df_subset)

    # B. Initialize Mapper (DB connection)
    mapper = GeneCoordinateMapper(get_paths()['db'])

    # C. Load raw genomic data
    gene_to_data = get_locus_to_data_dict(include_introns=True, gene_subset=genes_subset)

    # D. Build Sequence Registry (Pre-calculate mRNA sequences)
    registry = build_gene_sequence_registry(genes_subset, gene_to_data, mapper)

    # 3. EXECUTION
    # Note: df_subset must have 'sense_start' calculated before this step
    # (assuming it came from preprocess_aso_data or a previous step)
    df_result = add_external_mrna_and_context_columns(
        df=df_subset,
        mapper=mapper,
        gene_registry=registry
    )

    # 4. Clean & Check
    df_result = df_result.reindex(sorted(df_result.columns), axis=1)
    data_dict = df_result.to_dict(orient="records")
    data_regression.check(data_dict)