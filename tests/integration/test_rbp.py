import pytest

from tauso.features.rbp.load_rbp import load_attract_data
from tauso.populate.populate_rbp import encode_sequences, process_rbp


def test_per_motif_occupancy_regression(data_regression):
    """
    Comprehensive regression test for the per-motif scorer (the probability of occupying at
    least one site, 1 - exp(Sum log(1 - o))) using ALL real TAUSO PWM matrices against a highly
    varied sequence set (including overlapping motifs).
    """
    # 1. Load the real matrices
    rbp_map, pwm_db = load_attract_data()

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
    }

    # 2. Lay the sequences out the way the pipeline scans them
    names = list(test_sequences)
    flat_seq, offsets = encode_sequences(test_sequences.values())

    # 3. Collect results from the scorer
    results = {}

    # Sort the RBPs alphabetically to guarantee the YAML is strictly deterministic
    for rbp in sorted(rbp_map.keys()):
        matrix_ids = rbp_map[rbp]
        valid_mids = [m for m in matrix_ids if m in pwm_db]

        if not valid_mids:
            continue

        # Regression-test the scorer on one matrix; production noisy-OR-combines it over all of a protein's matrices
        _, occupancy = process_rbp({"matrices": [pwm_db[valid_mids[0]]], "col_name": rbp}, flat_seq, offsets)

        # Round to 8 decimal places to avoid cross-platform floating point drift
        results[rbp] = {name: round(float(value), 8) for name, value in zip(names, occupancy)}

    # 4. Check against the baseline
    # On the first run, this generates the baseline YAML file.
    # On future runs, it compares against it.
    data_regression.check(results)


def test_unknown_base_raises():
    """A sequence with an unknown base (N/etc.) must fail loudly, not be silently scored."""
    with pytest.raises(ValueError):
        encode_sequences(["AUGNNNUGA"])
