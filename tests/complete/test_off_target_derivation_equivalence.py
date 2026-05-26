"""
Validate the two assumptions behind collapsing off-target's 24 redundant RIsearch
passes (2 methods x 4 top_n x 3 cutoffs) into ONE loosest pass + in-memory
derivation:

  1. cutoff: a RIsearch run at a loose ``-s`` cutoff, filtered to ``score >
     strict_cutoff``, reproduces a run made directly at ``strict_cutoff``. Two
     load-bearing details, both verified here on real data: the cutoff filters
     on the SCORE column (not energy), and RIsearch's ``-s S`` is EXCLUSIVE
     (reports ``score > S``) — so the in-memory filter must be ``>``, not ``>=``
     (a hit at exactly ``score == S`` is absent from the direct ``-s S`` run).
  2. top_n: a RIsearch run against the top-N target genes, restricted to a
     top-K gene prefix (K < N), reproduces a run made directly against just the
     top-K genes (RIsearch scores each (query, target) pair independently, so
     dropping genes from the target FASTA cannot change the surviving rows).

Real genes (the ``gene_to_data`` session fixture), the real RIsearch binary, no
mocks. Queries are the antisense of real 30-mers drawn from the target mRNAs, so
each query has at least one strong self-hit that survives the strict cutoff —
the equivalence assertions are non-trivial.
"""

import os

import pandas as pd
import pytest

from tauso.features.hybridization.fast_hybridization import (
    TMP_PATH,
    Interaction,
    dump_target_file,
    risearch_hits_dataframe,
)
from tauso.util import get_antisense

LOOSE_CUTOFF = 800
STRICT_CUTOFF = 1200
SLICE_LEN = 4000
_SORT_COLS = ["trigger", "target", "trigger_start", "target_start", "score", "energy"]


def _run(query_pairs, target_path, cutoff):
    """One real RIsearch batch call -> raw hit DataFrame (8 columns)."""
    return risearch_hits_dataframe(
        query_pairs,
        target_path,
        minimum_score=cutoff,
        parsing_type="2",
        interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
        transpose=True,
    )


def _canonical(df):
    """Sort rows so two hit tables compare equal regardless of RIsearch output order."""
    return df[_SORT_COLS].sort_values(_SORT_COLS).reset_index(drop=True)


@pytest.fixture(scope="module")
def real_target_and_queries(gene_to_data):
    """3 real genes (mRNA sliced for speed) + one antisense 30-mer query per gene."""
    genes, seq_map = [], {}
    for gene, data in gene_to_data.items():
        mrna = data.full_mrna
        if mrna and len(mrna) >= SLICE_LEN:
            seq_map[gene] = mrna[:SLICE_LEN]
            genes.append(gene)
        if len(genes) == 3:
            break
    assert len(genes) == 3, "need 3 real genes with a long-enough mRNA"

    query_pairs = [(str(i), get_antisense(seq_map[g][500:530])) for i, g in enumerate(genes)]
    return genes, seq_map, query_pairs


def test_cutoff_filter_on_score_equals_direct_strict_run(real_target_and_queries):
    """Loose run filtered to score > strict == a run made directly at the strict cutoff."""
    genes, seq_map, query_pairs = real_target_and_queries
    TMP_PATH.mkdir(parents=True, exist_ok=True)
    target_path = dump_target_file("ot-deriv-cutoff.fa", seq_map)
    try:
        loose = _run(query_pairs, target_path, LOOSE_CUTOFF)
        strict = _run(query_pairs, target_path, STRICT_CUTOFF)
    finally:
        if os.path.exists(target_path):
            os.remove(target_path)

    assert not strict.empty, "strict run produced no hits — test would be vacuous"
    # The loose run must be a complete superset of the strict run (no strict hit
    # is missing from loose — RIsearch seeding does not drop hits as -s loosens).
    missing = (
        strict[_SORT_COLS]
        .merge(loose[_SORT_COLS], on=_SORT_COLS, how="left", indicator=True)
        .query("_merge == 'left_only'")
    )
    assert missing.empty, "loose run is not a superset of the strict run"
    # `-s` is exclusive, so the boundary is `>` (a hit at score == strict exists
    # in the loose run but not the strict run).
    derived = loose[loose["score"] > STRICT_CUTOFF]
    pd.testing.assert_frame_equal(_canonical(derived), _canonical(strict))


def test_topn_gene_prefix_filter_equals_direct_subset_run(real_target_and_queries):
    """Top-N run restricted to a top-K gene prefix == a run made directly against top-K."""
    genes, seq_map, query_pairs = real_target_and_queries
    prefix = genes[:2]
    TMP_PATH.mkdir(parents=True, exist_ok=True)
    full_path = dump_target_file("ot-deriv-full.fa", seq_map)
    sub_path = dump_target_file("ot-deriv-sub.fa", {g: seq_map[g] for g in prefix})
    try:
        full = _run(query_pairs, full_path, LOOSE_CUTOFF)
        sub = _run(query_pairs, sub_path, LOOSE_CUTOFF)
    finally:
        for p in (full_path, sub_path):
            if os.path.exists(p):
                os.remove(p)

    assert not sub.empty, "subset run produced no hits — test would be vacuous"
    derived = full[full["target"].isin(set(prefix))]
    pd.testing.assert_frame_equal(_canonical(derived), _canonical(sub))


def test_combined_derivation_equals_direct_strict_subset_run(real_target_and_queries):
    """The actual refactor path: ONE loose + all-genes run, filtered by both
    score > strict AND the top-K gene prefix, reproduces a direct strict + top-K run."""
    genes, seq_map, query_pairs = real_target_and_queries
    prefix = genes[:2]
    TMP_PATH.mkdir(parents=True, exist_ok=True)
    loose_full_path = dump_target_file("ot-deriv-comb-full.fa", seq_map)
    strict_sub_path = dump_target_file("ot-deriv-comb-sub.fa", {g: seq_map[g] for g in prefix})
    try:
        loose_full = _run(query_pairs, loose_full_path, LOOSE_CUTOFF)
        direct = _run(query_pairs, strict_sub_path, STRICT_CUTOFF)
    finally:
        for p in (loose_full_path, strict_sub_path):
            if os.path.exists(p):
                os.remove(p)

    assert not direct.empty, "direct run produced no hits — test would be vacuous"
    derived = loose_full[(loose_full["score"] > STRICT_CUTOFF) & (loose_full["target"].isin(set(prefix)))]
    pd.testing.assert_frame_equal(_canonical(derived), _canonical(direct))
