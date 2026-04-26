import os
import pytest
import pandas as pd

from notebooks.data.OligoAI.assign_canonical_gene import process_and_assign_genes

def generate_200_test_sequences():
    # Real transcript segments for high-confidence mapping
    seeds = {
        "KRAS": "AUGACUGAAUAUAAACUUGUGGUAGUUGGAGCUGGUGGCGUAGGCAAGAGUGCCUUGACGAUACAGCUAAUUCAGAAUCAUUUUGUGUAC",
        "SMN": "AUUGGUGGAGUGCUCCUCCUCCCCACCUCCCCUAUAUGUGUAGUCCCAGUUUUUUUUUUUUUUUUCUUUUAAAAGUAAAAGUCAUUUAG",
        "APOL1": "GAAUUAACUCCUACAGCAAUGGAGGGAGCUGCUCAUCUGAGAGCAGUCCUCCUCCUCCUCCUCCUCCUCCUCCUCCUCCUCCUCCUCCUC",
        "MYO6": "AUGGAGGACGAGAAGCUGAAGUUCGAGUUCGACGAGGAGUUCGAGCUCGAGUUCGACGAGGAGUUCGAGCUCGAGUUCGACGAGGAGUUC",
        "HTT": "AUGGCUAGUGCUAGUGCUAGUGCUAGUGCUAGUGCUAGUGCUAGUGCUAGUGCUAGUGCUAGUGCUAGUGCUAGUGCUAGUGCUAGUGCUA",
    }

    test_data = []

    # Generate 40 tiled sequences for each of the 5 genes (200 total)
    for gene, transcript in seeds.items():
        for i in range(40):
            # Extract a 20bp window
            seq = transcript[i : i + 20]
            if len(seq) == 20:
                test_data.append({"target_gene": gene, "rna_sequence": seq})

    # Add 10 "Scrambled" negatives at the end for controls
    negatives = [
        "GCACUACCAGAGCUAACUCA",  # eGFP
        "UGCAUACUUAACGCGUAUCG",  # Firefly Luciferase
        "CCGGAAUUUCGACGUGAAUA",  # LacZ
        "UACGUCGAUCAGUCGACUGC",  # Synthetic repeat
        "CGAAUACGGAUUACGCGAUG",  # Synthetic repeat
        "AUCGACGCUAGCGUAAUCGG",  # Synthetic repeat
    ]
    for neg in negatives:
        test_data.append({"target_gene": "UNKNOWN", "rna_sequence": neg})

    return test_data


TEST_DATA = generate_200_test_sequences()

def test_process_and_assign_genes_bulk_real_genome(tmp_path):
    """
    A true integration test hitting the real GRCh38 genome and Bowtie binary.
    No mocks are used.
    """

    # 1. Setup isolated file paths using pytest's tmp_path
    input_csv = tmp_path / "test_oligos_input.csv"
    output_csv = tmp_path / "test_oligos_output.csv"
    stats_csv = tmp_path / "test_oligos_output_STATS.csv"

    # 2. Create the physical input CSV
    df_input = pd.DataFrame(TEST_DATA)
    df_input.to_csv(input_csv, index=False)

    # 3. Execute the actual function
    # Note: Ensure genome="GRCh38" is installed in your env before running
    df_out, ambiguous_genes, df_stats = process_and_assign_genes(
        input_csv=str(input_csv),
        output_csv=str(output_csv),
        genome="GRCh38",
        sequence_col="rna_sequence",
    )

    # --- 4. ASSERTIONS ---

    # Check that the files were actually physically created
    assert os.path.exists(output_csv), "Output CSV was not created."
    assert os.path.exists(stats_csv), "Stats CSV was not created."

    # Check that the temp FASTA file was successfully cleaned up
    # (Your code uses a dynamic UUID, so we just ensure no fasta files are left in the test root)
    fasta_files = list(tmp_path.glob("*.fasta"))
    assert len(fasta_files) == 0, f"Temporary FASTA files were not deleted: {fasta_files}"

    # --- Validate Output DataFrames ---

    # Verify Stats DataFrame structure
    assert not df_stats.empty
    expected_stat_columns = [
        "original_target_gene",
        "total_asos_tested",
        "most_popular_canonical",
        "popular_hit_count",
        "match_rate",
        "all_hit_frequencies",
    ]
    for col in expected_stat_columns:
        assert col in df_stats.columns, f"Missing expected column in stats: {col}"

    # Verify the logic accurately counted our inputs
    # We fed it 40 KRAS sequences, 40 SMN, 40 APOL1, 40 MYO6, 6 UNKNOWN

    # Extract specific rows for testing
    kras_stats = df_stats[df_stats["original_target_gene"] == "KRAS"].iloc[0]
    smn_stats = df_stats[df_stats["original_target_gene"] == "SMN"].iloc[0]

    # Assert counts
    assert kras_stats["total_asos_tested"] == 40
    assert smn_stats["total_asos_tested"] == 40

    # Assuming Bowtie actually finds these in GRCh38, verify the canonical assignment worked
    assert kras_stats["most_popular_canonical"] == "KRAS"
    assert "40/40" in kras_stats["match_rate"]

    # Check the negative control (UNKNOWN gene)
    unknown_stats = df_stats[df_stats["original_target_gene"] == "UNKNOWN"].iloc[0]
    assert unknown_stats["total_asos_tested"] == 6
    assert pd.isna(unknown_stats["most_popular_canonical"]) # Should be NaN since it won't hit anything
    assert unknown_stats["popular_hit_count"] == 0
