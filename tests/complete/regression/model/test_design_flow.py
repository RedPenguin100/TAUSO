"""Full-circle integration test for the ASO design flow + brief input-validation tests.

The full-circle test tiles candidate ASOs across a target, featurizes them through the real pipeline,
scores them with the bundled efficacy model, and regression-checks the consumer view. Needs the full
feature-pipeline assets (genome, off-target rRNA, expression, ...); provisioned in CI by the setup-*
steps and locally by `tauso setup-*`. If an asset is missing it fails loudly rather than skipping.
The validation tests are fast (they raise before any featurization).
"""

import numpy as np
import pytest

from tauso.aso_generation import design_asos, design_details, summarize_design
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

    # consumer + safety views
    summary = summarize_design(ranked)
    assert list(summary["rank"]) == [1, 2, 3, 4, 5]
    details = design_details(ranked)
    assert {"flag_immune_cpg", "flag_binds_rrna", "flag_hepatotox_g4_grun", "liabilities"} <= set(details.columns)

    # regression-check the slim consumer view (stable columns; tolerance on the float score)
    dataframe_regression.check(summary, default_tolerance={"atol": 1e-4, "rtol": 1e-4})


@pytest.mark.integration
def test_design_asos_off_targets_table():
    """off_targets=True returns a tidy per-hit sequence off-target table next to the ranking.

    Needs the Bowtie index (`tauso setup-bowtie`), so it runs under --run-integration. TMSB10 is a
    custom gene here, so its own on-target gene symbol is passed via off_target_exclude_genes.
    """
    ranked, offtargets = design_asos(
        "USER_TMSB10",
        gene_sequence=TMSB10_PREMRNA,
        first_n=5,
        top_n=5,
        off_targets=True,
        off_target_max_distance=2,
        off_target_exclude_genes=["TMSB10"],
        n_jobs=1,
    )

    assert len(ranked) == 5
    assert list(offtargets.columns) == [
        "rank",
        "aso_sequence",
        "off_target_gene",
        "distance",
        "region",
        "chrom",
        "start",
        "strand",
    ]
    if not offtargets.empty:
        assert offtargets["distance"].between(0, 2).all()  # exact distance within the requested range
        assert set(offtargets["aso_sequence"]).issubset(set(ranked[ASO_SEQUENCE]))
        assert offtargets["rank"].between(1, 5).all()
        assert "TMSB10" not in set(offtargets["off_target_gene"])  # on-target gene excluded


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


# The ten top-ranked candidates from the full MALAT1 design (A549 / Lipofection / 100 nM, 20-mer
# gapmer), in descending predicted-efficacy order. Tiling + scoring all ~8,800 MALAT1 candidates takes
# ~24 min, but every feature is computed per-candidate against the gene (not relative to the scored
# batch), so re-scoring exactly these ten in isolation reproduces the full-run scores bit-for-bit.
# Freezing them guards the published ranking with a fast, always-on test and a 10-row fixture.
MALAT1_TOP10 = (
    "GCCACTTCCTTTGCTCTGCA",
    "TCTCATTTATTTCGGCTTCT",
    "CCTTAGTTGGCATCAAGGCA",
    "TCCTCATTTTTTCAGTGCTA",
    "TCTCTCATTTATTTCGGCTT",
    "CTCTCATTTATTTCGGCTTC",
    "GCTATTTTTTCTTACTGGGT",
    "ATCTCTCATTTATTTCGGCT",
    "GGGTTATGCTTATTCCCCAA",
    "CGGCTTCTTTTATTCCAGGA",
)


def test_malat1_top10_frozen(dataframe_regression):
    """Freeze the published MALAT1 top-10 ASO ranking against model/feature/data drift.

    Re-scores exactly the ten frozen candidates through the real feature + model pipeline and
    regression-checks the consumer view. This is the same computation `design_asos` runs on the full
    gene, restricted to these ten candidates (a valid shortcut: features are per-candidate), so the
    scores and rank order must reproduce the full design.
    """
    from tauso.aso_generation import (
        AssetCache,
        Transfection,
        _apply_standard_metadata,
        _fill_out_of_range_one_hots,
        default_config,
        generate_aso_features,
        get_initial_data,
    )
    from tauso.inference import load_model, predict

    config = default_config()
    config.cell_line = "A549"
    config.transfection_method = Transfection.LIPOFECTION
    config.volume = 100
    config.organism_name = "human"

    cache = AssetCache(genome="GRCh38")
    malat1 = cache.get_full_gene_data()["MALAT1"].full_mrna
    candidates = get_initial_data(malat1, aso_sizes=[20], canonical_name="MALAT1")
    frozen = candidates[candidates[ASO_SEQUENCE].isin(MALAT1_TOP10)].copy()
    assert len(frozen) == len(MALAT1_TOP10), "a frozen top-10 sequence is no longer a MALAT1 20-mer candidate"

    _apply_standard_metadata(frozen, config)
    featured, _ = generate_aso_features(frozen, cache, n_jobs=1)
    _, model_features = load_model()
    _fill_out_of_range_one_hots(featured, model_features)

    col = score_column()
    featured[col] = predict(featured)
    ranked = featured.sort_values(col, ascending=False, kind="stable").reset_index(drop=True)

    summary = summarize_design(ranked)
    assert list(summary["rank"]) == list(range(1, len(MALAT1_TOP10) + 1))
    assert list(summary[ASO_SEQUENCE]) == list(MALAT1_TOP10)  # score order matches the published ranking
    dataframe_regression.check(summary, default_tolerance={"atol": 1e-4, "rtol": 1e-4})
