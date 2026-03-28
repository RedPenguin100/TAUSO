import os

import pandas as pd
import pytest
from pandarallel import pandarallel

from tauso.data.consts import SEQUENCE
from tauso.features.hybridization_off_target.off_target_specific_gene import (
    off_target_specific_seq_pandarallel,
)
from tests.conftest import SHORT_GENE

# Initialize pandarallel (run this once, e.g., at the top of your notebook/script)
# progress_bar=True is very helpful for long jobs
pandarallel.initialize(progress_bar=True, nb_workers=os.cpu_count())


@pytest.mark.integration
def test_off_target_regression(
    small_gene_to_data, dataframe_regression, data_regression
):
    """
    Regression test for off_target_specific_seq.
    Checks if the output DF and feature name match the stored snapshot.
    """
    sequences = [
        # --- Original Set (Reference 20-mers) ---
        "ATCGATCGATCGATCGATCG",  # 1. Simple repeating
        "GGGGCCCCTTTTAAAAATCG",  # 2. Block pattern
        "TGCATGCATGCATGCATGCA",  # 3. Mixed content
        "AAAAAAAAAAAAAAAAAAAA",  # 4. Homopolymer (Poly-A, 20-mer)
        # --- Homopolymers & Simple Repeats (Edge Cases) ---
        "TTTTTTTTTTTTTTT",  # 5. Poly-T (15-mer)
        "CCCCCCCCCCCCCCCCCC",  # 6. Poly-C (18-mer)
        "GGGGGGGGGGGGGGGGGGGGGG",  # 7. Poly-G (22-mer)
        "ATATATATATATAT",  # 8. Dinucleotide repeat AT (14-mer)
        "GCGCGCGCGCGCGCGCGC",  # 9. Dinucleotide repeat GC (18-mer)
        "CAGCAGCAGCAGCAGCAG",  # 10. Trinucleotide repeat CAG (18-mer)
        "TGTGTGTGTGTGTGTGTGTG",  # 11. TG repeat (20-mer)
        # --- Low GC Content (AT-Rich) ---
        "AAATTTATATTTTAA",  # 12. Low GC (15-mer)
        "TATATAAATTTTAAAT",  # 13. Low GC (16-mer)
        "ATTTTAAAATTTTAAAATTT",  # 14. Low GC (20-mer)
        "AAATTTCCAAATTT",  # 15. Mostly AT with sparse C (14-mer)
        "TTTAAATTTGGTTTAAA",  # 16. Mostly AT with sparse G (17-mer)
        "ATAATAATAATAATAATAAT",  # 17. AT only (20-mer)
        "TTAATTAATTAATTAA",  # 18. TTAA repeats (16-mer)
        # --- High GC Content (GC-Rich) ---
        "GGCCGGCCGGCCGG",  # 19. High GC (14-mer)
        "GCGCGGCCGGCGCGGCC",  # 20. High GC (17-mer)
        "CCGGCCGGCCGGCCGGCCGG",  # 21. High GC (20-mer)
        "GGGGCCCCGGGGCCCCGGGG",  # 22. High GC (20-mer)
        "CGCGCGCGCGTGCGCGCGCG",  # 23. High GC with single T (20-mer)
        "GGCCAAAAGGCC",  # 24. High GC short (12-mer? -> adjusted to 14) -> "GGCCAAAAGGCCGG"
        "CGGCGGCCGGCGGCGGCCGGCC",  # 25. Very High GC (22-mer)
        # --- Varying Lengths (Random Mixed) ---
        "AGCTAGCTAGCTAG",  # 26. 14-mer
        "TCGATCGATCGATCG",  # 27. 15-mer
        "GCTAGCTAGCTAGCTA",  # 28. 16-mer
        "CTAGCTAGCTAGCTAGC",  # 29. 17-mer
        "TAGCTAGCTAGCTAGCTA",  # 30. 18-mer
        "AGCTAGCTAGCTAGCTAGC",  # 31. 19-mer
        "GATCGATCGATCGATCGATC",  # 32. 20-mer
        "ATCGATCGATCGATCGATCGA",  # 33. 21-mer
        "TCGATCGATCGATCGATCGATC",  # 34. 22-mer
        # --- Structural Motifs ---
        "GAATTCGAATTCGAATTC",  # 35. EcoRI site repeats (18-mer)
        "GGATCCGGATCCGGATCC",  # 36. BamHI site repeats (18-mer)
        "TTTTTGGGGGAAAAACCCCC",  # 37. Blocks of 5 (20-mer)
        "ACGTACGTACGTAC",  # 38. Tetranucleotide repeat (14-mer)
        # --- Random Biological-Like Sequences ---
        "CCAAGGTTCCGGTTAA",  # 39. Mixed (16-mer)
        "AACCGGTTTTGGCCAA",  # 40. Mixed (16-mer)
        "TGACGTACGTTAGC",  # 41. Mixed (14-mer)
        "CTGACGTACGTTTAGCT",  # 42. Mixed (17-mer)
        "ACTGACGTACGTTTAGCTG",  # 43. Mixed (19-mer)
        "GACTGACGTACGTTTAGCTGA",  # 44. Mixed (21-mer)
        "TGACTGACGTACGTTTAGCT",  # 45. Mixed (20-mer)
        "GTACGTTTAGCTGACGTACG",  # 46. Mixed (20-mer)
        "TACGTTTAGCTGACGTACGT",  # 47. Mixed (20-mer)
        "ACGTTTAGCTGACGTACGTT",  # 48. Mixed (20-mer)
        "CGTTTAGCTGACGTACGTTT",  # 49. Mixed (20-mer)
        "GTTTAGCTGACGTACGTTTA",  # 50. Mixed (20-mer)
    ]

    # Create the DataFrame with actual data
    aso_df = pd.DataFrame({SEQUENCE: sequences})
    aso_df = aso_df.reset_index()

    # 2. Execution
    result_df, feat_name = off_target_specific_seq_pandarallel(
        aso_df, SHORT_GENE, small_gene_to_data, cutoff=800, n_jobs=1
    )

    # 3. Regression Checks

    # Check the main DataFrame.
    # This will create a CSV file on the first run, and compare against it on future runs.
    dataframe_regression.check(result_df, basename="off_target_df_regression")

    # Check the feature name (scalar string).
    # This ensures the secondary return value hasn't changed.
    data_regression.check(
        {"feature_name": feat_name}, basename="off_target_meta_regression"
    )
