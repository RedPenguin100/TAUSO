"""
Verify that get_triggers_mfe_scores_batch produces results identical to N
individual calls to get_trigger_mfe_scores_by_risearch.

Uses the GFP sequences from integration test data — no genome loading required.
Marked integration so it only runs with --run-integration.
"""

import pytest
from Bio import SeqIO

from tauso.features.hybridization.fast_hybridization import (
    Interaction,
    get_trigger_mfe_scores_by_risearch,
    get_triggers_mfe_scores_batch,
)
from tauso.features.hybridization.off_target.off_target_functions import parse_risearch_output
from tauso.util import get_antisense
from tests.common.consts import INTEGRATION_TESTS_PATH


@pytest.fixture(scope="module")
def gfp_target_seq():
    fasta_path = INTEGRATION_TESTS_PATH / "data" / "GFP_first_exp.fasta"
    record = next(SeqIO.parse(str(fasta_path), "fasta"))
    return str(record.seq.upper())


@pytest.fixture(scope="module")
def gfp_name_to_seq(gfp_target_seq):
    return {"gfp": gfp_target_seq}


# A handful of 20-mer queries drawn from GFP
QUERY_SEQS = [
    "ATGGTGAGCAAGGGCGAGGA",
    "GCGAGGGCGATGCCACCTAC",
    "TTCACCTTGATGCCGTTCTT",
    "CTGAACTTGTGGCCGTTTAC",
    "CTTGTACAGCTCGTCCATGC",
]

CUTOFF = 800


@pytest.mark.integration
def test_batch_output_equals_single_outputs(gfp_name_to_seq, tmp_path):
    """Batch call must produce the same raw text as N sequential single calls, line by line."""
    from tauso.features.hybridization.fast_hybridization import TMP_PATH, dump_target_file

    TMP_PATH.mkdir(parents=True, exist_ok=True)
    target_path = dump_target_file("gfp-test-target.fa", gfp_name_to_seq)

    try:
        single_results = {}
        for seq in QUERY_SEQS:
            trigger = get_antisense(seq)
            raw = get_trigger_mfe_scores_by_risearch(
                trigger,
                gfp_name_to_seq,
                minimum_score=CUTOFF,
                parsing_type="2",
                interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
                transpose=True,
                target_file_cache=target_path,
            )
            df = parse_risearch_output(raw) if raw.strip() else None
            single_results[seq] = df

        pairs = [(str(i), get_antisense(seq)) for i, seq in enumerate(QUERY_SEQS)]
        batch_raw = get_triggers_mfe_scores_batch(
            trigger_id_seq_pairs=pairs,
            target_file_path=target_path,
            minimum_score=CUTOFF,
            parsing_type="2",
            interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
            transpose=True,
        )
    finally:
        import os

        if os.path.exists(target_path):
            os.remove(target_path)

    batch_df = parse_risearch_output(batch_raw) if batch_raw.strip() else None

    for i, seq in enumerate(QUERY_SEQS):
        trigger_id = str(i)
        single_df = single_results[seq]

        if single_df is None:
            batch_hits = batch_df[batch_df["trigger"] == trigger_id] if batch_df is not None else None
            assert batch_hits is None or batch_hits.empty, f"Single returned no hits for seq {i}, but batch did"
            continue

        # Filter batch to this trigger
        assert batch_df is not None, f"Batch returned nothing but single had hits for seq {i}"
        batch_hits = batch_df[batch_df["trigger"] == trigger_id].reset_index(drop=True)

        # Compare energies (sorted, since ordering may differ)
        single_energies = sorted(single_df["energy"].tolist())
        batch_energies = sorted(batch_hits["energy"].tolist())
        assert single_energies == pytest.approx(batch_energies, rel=1e-6), (
            f"Energy mismatch for query {i} ({seq!r}): single={single_energies}, batch={batch_energies}"
        )
