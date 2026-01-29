import pytest
import pandas as pd

# Imports from your project
from notebooks.consts import UPDATED_CSV
from notebooks.preprocessing import preprocess_aso_data, get_unique_genes
from tauso.data.data import get_paths
from tauso.genome.TranscriptMapper import GeneCoordinateMapper, build_gene_sequence_registry
from tauso.genome.read_human_genome import get_locus_to_data_dict
from tauso.algorithms.genomic_context_windows import add_external_mrna_and_context_columns

from tauso.data.data import get_data_dir, load_db
from tauso.features.codon_usage.find_cai_reference import filter_gtf_genes


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


def test_real_genome_filter():
    """
    Connects to the REAL local GRCh38 database and verifies
    that specific biological cases are filtered correctly.
    """
    # 1. Load Real DB (Change genome version if you use hg19)
    try:
        # Assuming load_db is available in your environment
        db = load_db(genome="GRCh38")
    except Exception as e:
        pytest.fail(f"Could not load real database: {e}")

    print(f"\nConnected to DB: {db.dbfn}")

    # --- Case 1: Strict Protein Coding ---
    print("Testing 'protein_coding' mode...")
    strict_genes = filter_gtf_genes(db, filter_mode='protein_coding')

    # 1. Standard Protein Coding -> MUST BE IN
    assert "GAPDH" in strict_genes, "GAPDH (Protein Coding) is missing!"
    assert "TP53" in strict_genes, "TP53 (Protein Coding) is missing!"

    # 2. lncRNA -> MUST BE OUT
    assert "MALAT1" not in strict_genes, "Error: MALAT1 (lncRNA) was NOT filtered out!"
    assert "NEAT1" not in strict_genes, "Error: NEAT1 (lncRNA) was NOT filtered out!"

    # 3. Mitochondrial -> MUST BE OUT
    # Note: Check your DB to see if it uses 'MT-ND1' or 'ND1'
    assert "MT-ND1" not in strict_genes, "Error: MT-ND1 (Mitochondrial) was NOT filtered out!"

    # --- Case 2: Loose Non-Mitochondrial ---
    print("Testing 'non_mt' mode...")
    loose_genes = filter_gtf_genes(db, filter_mode='non_mt')

    # 1. lncRNA -> MUST BE IN
    assert "MALAT1" in loose_genes, "Error: MALAT1 should be present in non_mt mode!"

    # 2. Mitochondrial -> MUST BE OUT
    assert "MT-ND1" not in loose_genes, "Error: MT-ND1 should still be filtered in non_mt mode!"

    print("\n✅ Real DB Filter Test Passed!")