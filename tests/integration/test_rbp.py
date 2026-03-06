import numpy as np
import pandas as pd
import pytest
from numba import njit

# Import your actual data loaders and the original function
from tauso.features.rbp.load_rbp import load_attract_data
from tauso.features.rbp.pwm_helper import calculate_total_affinity

import numpy as np
import pandas as pd
import pytest

# Import your actual data loaders and the affinity function
from tauso.features.rbp.load_rbp import load_attract_data
from tauso.populate.populate_rbp import calculate_total_affinity_numba


# Note: Import your original function FIRST to generate the baseline YAML,
# then swap to calculate_total_affinity_numba to verify it matches.


def test_calculate_total_affinity_regression(data_regression):
    """
    Comprehensive regression test for affinity math using ALL real TAUSO PWM matrices
    against a highly varied sequence set (including overlapping motifs).
    """
    # 1. Load the real matrices
    rbp_map, pwm_db = load_attract_data()

    # 2. Define 20 test sequences covering motifs, overlaps, repeats, and edge cases
    test_sequences = {
        # Overlapping & well-known biological motifs
        "fox_overlap": "UGCAUGCAUGCAUGCAUGCA",  # RBFOX overlapping motifs
        "pumilio_motif": "UGUAAUAUUGUAAUAUUGUA",  # PUF family motifs
        "srsf_rich": "GGAGGAGGAGGAGGAGGAGG",  # SRSF targets (Purine rich)
        "hnrnp_rich": "UAGUAGUAGUAGUAGUAGUA",  # hnRNP targets
        "cug_repeat": "CUGCUGCUGCUGCUGCUGCU",  # CELF/MBNL overlapping targets
        "au_rich_element": "UUAUUUAUUAUUUAUUAUUU",  # AREs (e.g., HuR/ELAVL1 binding)
        "gc_rich_hairpin": "GCGCGCGCUUUUGCGCGCGC",  # GC clamp with a U loop

        # Homopolymers (tests extreme affinity for specific RBPs like PolyA Binding Protein)
        "poly_A": "AAAAAAAAAAAAAAAAAAAA",
        "poly_U": "UUUUUUUUUUUUUUUUUUUU",
        "poly_C": "CCCCCCCCCCCCCCCCCCCC",
        "poly_G": "GGGGGGGGGGGGGGGGGGGG",

        # Complex / Mixed
        "mixed_complex_1": "AUGUCGACGUUAGCAUGCUA",
        "mixed_complex_2": "CGCGCGAUAUAUAGCGCGCG",
        "alternating_ry": "CUCUCUCUCUCUCUCUCUCU",  # Alternating Pyrimidine/Purine

        # Formatting / Edge Cases
        "dna_version": "TGCATGCATGCATGCATGCA",  # Ensures 'T' maps to 'U' properly
        "lowercase_seq": "ugcauguauuauggag",  # Ensures case insensitivity
        "short_seq": "AUG",  # Should trigger the 0.0 short circuit
        "empty_seq": "",  # Should trigger the 0.0 empty circuit
        "nan_seq": np.nan,  # Pandas missing data
        "dirty_seq": "AUGNNNUGAXXXC"  # Unknown bases should be ignored/skipped
    }

    # Custom background probs to ensure that logic path is tested (e.g., A/T rich background)
    custom_bg = np.array([0.3, 0.2, 0.2, 0.3])

    # 3. Collect results from the function
    results = {}

    # Sort the RBPs alphabetically to guarantee the YAML is strictly deterministic
    for rbp in sorted(rbp_map.keys()):
        matrix_ids = rbp_map[rbp]
        valid_mids = [m for m in matrix_ids if m in pwm_db]

        if not valid_mids:
            continue

        # Use the first valid matrix exactly like your main code does
        real_pwm_matrix = pwm_db[valid_mids[0]]

        # Sub-dictionary for this specific RBP
        rbp_results = {}

        for seq_name, seq_val in test_sequences.items():
            # Test with default background (None)
            score_default = calculate_total_affinity_numba(
                sequence=seq_val,
                pwm_matrix=real_pwm_matrix,
                background_probs=None
            )
            # Test with custom background
            score_custom = calculate_total_affinity_numba(
                sequence=seq_val,
                pwm_matrix=real_pwm_matrix,
                background_probs=custom_bg
            )

            # Round to 8 decimal places to avoid cross-platform floating point drift
            rbp_results[f"{seq_name}_default_bg"] = round(float(score_default), 8)
            rbp_results[f"{seq_name}_custom_bg"] = round(float(score_custom), 8)

        # Assign the sub-dictionary to the main results dictionary
        results[rbp] = rbp_results

    # 4. Check against the baseline
    # On the first run, this generates the baseline YAML file.
    # On future runs, it compares against it.
    data_regression.check(results)