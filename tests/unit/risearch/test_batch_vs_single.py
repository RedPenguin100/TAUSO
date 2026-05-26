"""
Verify a single batched RIsearch call yields the same per-query hits as N separate
single-query calls (the batching invariant: putting many queries in one FASTA must
not change any query's results).

Uses the GFP sequence from integration test data — no genome loading required.
Marked integration so it only runs with --run-integration.
"""

import os

import pytest
from Bio import SeqIO

from tauso.features.hybridization.fast_hybridization import (
    TMP_PATH,
    Interaction,
    dump_target_file,
    risearch_hits_dataframe,
)
from tauso.util import get_antisense
from tests.common.consts import INTEGRATION_TESTS_PATH


@pytest.fixture(scope="module")
def gfp_name_to_seq():
    fasta_path = INTEGRATION_TESTS_PATH / "data" / "GFP_first_exp.fasta"
    record = next(SeqIO.parse(str(fasta_path), "fasta"))
    return {"gfp": str(record.seq.upper())}


# A handful of 20-mer queries drawn from GFP
QUERY_SEQS = [
    "ATGGTGAGCAAGGGCGAGGA",
    "GCGAGGGCGATGCCACCTAC",
    "TTCACCTTGATGCCGTTCTT",
    "CTGAACTTGTGGCCGTTTAC",
    "CTTGTACAGCTCGTCCATGC",
]

CUTOFF = 800


def _hits(pairs, target_path):
    return risearch_hits_dataframe(
        pairs,
        target_path,
        minimum_score=CUTOFF,
        parsing_type="2",
        interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
        transpose=True,
    )


@pytest.mark.integration
def test_batch_hits_equal_single_query_hits(gfp_name_to_seq):
    """Batching N queries in one call must give each query the same hits as a solo call."""
    TMP_PATH.mkdir(parents=True, exist_ok=True)
    target_path = dump_target_file("gfp-test-target.fa", gfp_name_to_seq)
    try:
        pairs = [(str(i), get_antisense(seq)) for i, seq in enumerate(QUERY_SEQS)]
        batch = _hits(pairs, target_path)
        singles = {str(i): _hits([(str(i), get_antisense(seq))], target_path) for i, seq in enumerate(QUERY_SEQS)}
    finally:
        if os.path.exists(target_path):
            os.remove(target_path)

    for i in range(len(QUERY_SEQS)):
        tid = str(i)
        batch_energies = sorted(batch[batch["trigger"] == tid]["energy"].tolist())
        single_energies = sorted(singles[tid]["energy"].tolist())
        assert batch_energies == pytest.approx(single_energies, rel=1e-6), (
            f"Energy mismatch for query {i}: batch={batch_energies}, single={single_energies}"
        )
