import pytest
from Bio import SeqIO

from tauso.features.hybridization.fast_hybridization import (
    Interaction,
    get_mfe_scores,
    get_trigger_mfe_scores_by_risearch,
)
from tauso.util import get_antisense
from tests.common.consts import INTEGRATION_TESTS_PATH


def get_gfp_seq_and_context():
    gfp_context_path = INTEGRATION_TESTS_PATH / "data" / "GFP_context.txt"
    gfp_first_exp_path = INTEGRATION_TESTS_PATH / "data" / "GFP_first_exp.fasta"

    gfp_obj = next(SeqIO.parse(str(gfp_first_exp_path), "fasta"))
    gfp_seq = str(gfp_obj.seq.upper())

    with open(str(gfp_context_path), "r") as f:
        gfp_context = f.read().upper()

    gfp_start = gfp_context.find(gfp_seq)
    if gfp_start == -1:
        raise ValueError("Context not found!")

    return (gfp_seq, gfp_context)


def get_gfp_first_exp(gap=100):
    # TODO: gap should be always 100 in this function
    gfp_seq, gfp_context = get_gfp_seq_and_context()

    gfp_start = gfp_context.find(gfp_seq)
    if gfp_start == -1:
        raise ValueError("Context not found!")

    gfp_ext = gfp_context[gfp_start - gap : gfp_start + len(gfp_seq) + gap]

    return gfp_ext


def get_gfp_second_exp():
    right_gap = 50
    gfp_seq, gfp_context = get_gfp_seq_and_context()

    gfp_start = gfp_context.find(gfp_seq)
    gfp_ext = gfp_context[gfp_start : gfp_start + len(gfp_seq) + right_gap]

    return gfp_ext


@pytest.mark.integration
def test_risearch_gfp_modified_transpose(data_regression):
    # --- 2. Data Preparation ---
    gfp_seq = get_gfp_first_exp(gap=0)
    sample_seq = gfp_seq[:20]
    name_to_seq = {f"gfp_seq": gfp_seq}

    # Initialize a dictionary to hold all data we want to snapshot
    regression_data = {}

    # Capture context data (replacing your initial print statements)
    regression_data["metadata"] = {
        "gfp_ontarget_slice": gfp_seq[695:714],
        "sample_seq": sample_seq,
        "sample_antisense": get_antisense(sample_seq),
    }

    # --- 3. Test "Good" Sample ---
    result_good = get_trigger_mfe_scores_by_risearch(
        sample_seq,
        name_to_seq,
        interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
        minimum_score=900,
        neighborhood=30,
        parsing_type="2",
        transpose=True,
    )
    mfe_scores_good = get_mfe_scores(result_good, "2")

    print(mfe_scores_good)

    # Add results to our data snapshot
    regression_data["good_sample_run"] = {
        "risearch_result": result_good,
        "mfe_scores": mfe_scores_good,
    }

    # --- 4. Test "Bad" Samples ---
    bad_samples = [s + sample_seq[3:20] for s in ["AAA", "ATA", "AGA", "ACG"]]

    regression_data["bad_sample_runs"] = []

    for bad_sample in bad_samples:
        result_bad = get_trigger_mfe_scores_by_risearch(
            bad_sample,
            name_to_seq,
            interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
            minimum_score=900,
            neighborhood=30,
            parsing_type="2",
            transpose=True,
        )
        mfe_scores_bad = get_mfe_scores(result_bad, "2")

        # Append each iteration's result to the list
        regression_data["bad_sample_runs"].append(
            {
                "input_sequence": bad_sample,
                "risearch_result": result_bad,
                "mfe_scores": mfe_scores_bad,
            }
        )

    # --- 5. Perform Regression Check ---
    # This checks the current `regression_data` dict against the stored YAML file.
    data_regression.check(regression_data)


@pytest.mark.integration
def test_risearch_gfp_modified_original(data_regression):
    gfp_seq = get_gfp_first_exp(gap=0)
    sample_seq = gfp_seq[:20]
    name_to_seq = {f"gfp_seq": gfp_seq}

    # Initialize a dictionary to hold all data we want to snapshot
    regression_data = {}

    # Capture context data (replacing your initial print statements)
    regression_data["metadata"] = {
        "gfp_ontarget_slice": gfp_seq[695:714],
        "sample_seq": sample_seq,
        "sample_antisense": get_antisense(sample_seq),
    }

    # --- 3. Test "Good" Sample ---
    result_good = get_trigger_mfe_scores_by_risearch(
        sample_seq,
        name_to_seq,
        interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
        minimum_score=900,
        neighborhood=30,
        parsing_type="2",
        transpose=False,
    )
    mfe_scores_good = get_mfe_scores(result_good, "2")

    print(mfe_scores_good)

    # Add results to our data snapshot
    regression_data["good_sample_run"] = {
        "risearch_result": result_good,
        "mfe_scores": mfe_scores_good,
    }

    # --- 4. Test "Bad" Samples ---
    bad_samples = [s + sample_seq[3:20] for s in ["AAA", "ATA", "AGA", "ACG"]]

    regression_data["bad_sample_runs"] = []

    for bad_sample in bad_samples:
        result_bad = get_trigger_mfe_scores_by_risearch(
            bad_sample,
            name_to_seq,
            interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
            minimum_score=900,
            neighborhood=30,
            parsing_type="2",
            transpose=False,
        )
        mfe_scores_bad = get_mfe_scores(result_bad, "2")

        # Append each iteration's result to the list
        regression_data["bad_sample_runs"].append(
            {
                "input_sequence": bad_sample,
                "risearch_result": result_bad,
                "mfe_scores": mfe_scores_bad,
            }
        )

    # --- 5. Perform Regression Check ---
    # This checks the current `regression_data` dict against the stored YAML file.
    data_regression.check(regression_data)
