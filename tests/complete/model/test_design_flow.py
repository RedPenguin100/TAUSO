"""Full-circle integration test for the ASO design flow + brief input-validation tests.

The full-circle test tiles candidate ASOs across a target, featurizes them through the real pipeline,
scores them with the bundled efficacy model, and regression-checks the consumer view. Needs the full
feature-pipeline assets (genome, off-target rRNA, expression, ...); provisioned in CI by the setup-*
steps and locally by `tauso setup-*`. If an asset is missing it fails loudly rather than skipping.
The validation tests are fast (they raise before any featurization).
"""

import numpy as np
import pytest

from tauso.aso_generation import design_asos, feature_details, summarize_design, tox_details
from tauso.data.consts import ASO_SEQUENCE
from tauso.inference import score_column

TMSB10_PREMRNA = "GTTTCTTGCTGCAGCAACGCGAGTGGGAGCACCAGGATCTCGGGCTCGGAACGAGACTGCACGGTGAGTGCGGCGCCGGGGCGGGGGGCCCACCCAGGGTGTGGTCGGATCCGGTGCACCGGGCGGGCGCGCCGCAACCGCGACAGGCGCCCTTCTCGGACCGGACGCAGGGGCCGGCGACCACGCCCTGGGACCGAGAAGAGGGGTGCGGGACGCGCCCAGATCCTCGGCCTTGGGGCTGCTCGGCACGCCTTGGCGCGAGTGCCACGTCGAGAGGCGTCGGCGGGGAGCGCGGAAGGGGACGCTGGCCCCCAGGCCCAGGTCAAGCGCCTTGGTTTGCCCACTAGGATTGTTTTAAGAAAATGGCAGACAAACCAGACATGGGGGAAATCGCCAGCTTCGATAAGGCCAAGCTGAAGAAAACGGAGACGCAGGAGAAGAACACCCTGCCGACCAAAGAGAGTGAGTGTGCCTCGGTCTCCCGCGCCCCAGCCCAGCCCCTCACCCTGCTCTTCCTTGCAAACCCACTCCTCCACCCCCCACCCCGCCGTTGTCCCCGGTGTGGGCGGCCCCGGCCACTCTTTCAGTTTCACAAAGCGCCTTGTTTCTCCCCAGCCCCAAGCTTCCTTCTAAATCCCCACACCTCGTGGGTGCCTCGCCCACACCGGGAAGCACCTCGGTTGCGGGTGGGGGTTGCAGCTCCCCTCCAGCGCCCGCTTCCCGCTCTCCACAGCCATTGAGCAGGAGAAGCGGAGTGAAATTTCCTAAGATCCTGGAGGATTTCCTACCCCCGTCCTCTTCGAGACCCCAGTCGTGATGTGGAGGAAGAGCCACCTGCAAGATGGACACGAGCCACAAGCTGCACTGTGAACCTGGGCACTCCGCGCCGATGCCACCGGCCTGTGGGTCTCTGAAGGGACCCCCCCCCAATCGGACTGCCAAATTCTCCGGTTTGCCCCGGGATATTATAGAAAATTATTTGTATGAATAATGAAAATAAAACACACCTCGTGGCA"


def test_design_asos_full_circle(dataframe_regression):
    ranked = design_asos("USER_TMSB10", gene_sequence=TMSB10_PREMRNA, first_n=5, top_n=5, n_jobs=1)

    col = score_column()
    assert len(ranked) == 5
    assert col in ranked.columns
    scores = ranked[col].to_numpy()
    assert np.isfinite(scores).all()
    assert list(scores) == sorted(scores, reverse=True)  # ranked best-first
    assert ranked[ASO_SEQUENCE].str.len().eq(20).all()  # 20-mer candidates

    # consumer + tox + feature views (each keyed by the chemistry columns + sequence)
    chem = {"chemical_pattern", "modification", "ps_pattern"}
    summary = summarize_design(ranked)
    assert list(summary["rank"]) == [1, 2, 3, 4, 5]
    assert chem <= set(summary.columns)
    tox = tox_details(ranked)
    assert (chem | {"flag_immune_cpg", "flag_hepatotox_g4_grun", "liabilities"}) <= set(tox.columns)
    feats = feature_details(ranked)
    assert (chem | {"offtarget_transcriptome", "offtarget_rrna", "rnaseh1_cleavage_fit"}) <= set(feats.columns)

    # regression-check the slim consumer view (stable columns; tolerance on the float score)
    dataframe_regression.check(summary, default_tolerance={"atol": 1e-4, "rtol": 1e-4})


@pytest.mark.parametrize(
    "kwargs, match",
    [
        (dict(gene_sequence="ACGT", aso_sizes=[]), "non-empty"),
        (dict(gene_sequence="ACGT", aso_sizes=[30]), "within"),  # over the 28 nt upper bound
        (dict(gene_sequence="ACGT", aso_sizes=[10]), "within"),  # under the 12 nt lower bound
        (dict(gene_sequence="ACG ACG"), "whitespace"),
        (dict(gene_sequence="ACGZT"), "non-ACGU"),
        (dict(gene_sequence="ACGTACGTACGTACGT", aso_sizes=[20]), "exceeds target length"),
    ],
)
def test_design_asos_input_validation(kwargs, match):
    with pytest.raises(ValueError, match=match):
        design_asos("USER_TEST", **kwargs)
