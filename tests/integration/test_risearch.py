import pytest

from tauso.hybridization.fast_hybridization import get_trigger_mfe_scores_by_risearch, Interaction, get_mfe_scores
from tauso.util import get_antisense
from tests.integration.target_finder import get_gfp_first_exp


def test_risearch_gfp_modified_transpose(data_regression):
    # --- 2. Data Preparation ---
    gfp_seq = get_gfp_first_exp(gap=0)
    sample_seq = gfp_seq[:20]
    name_to_seq = {f"gfp_seq": gfp_seq}

    # Initialize a dictionary to hold all data we want to snapshot
    regression_data = {}

    # Capture context data (replacing your initial print statements)
    regression_data['metadata'] = {
        'gfp_ontarget_slice': gfp_seq[695:714],
        'sample_seq': sample_seq,
        'sample_antisense': get_antisense(sample_seq)
    }

    # --- 3. Test "Good" Sample ---
    result_good = get_trigger_mfe_scores_by_risearch(
        sample_seq,
        name_to_seq,
        interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
        minimum_score=900,
        neighborhood=30,
        parsing_type='2',
        transpose=True
    )
    mfe_scores_good = get_mfe_scores(result_good, '2')

    print(mfe_scores_good)

    # Add results to our data snapshot
    regression_data['good_sample_run'] = {
        'risearch_result': result_good,
        'mfe_scores': mfe_scores_good
    }

    # --- 4. Test "Bad" Samples ---
    bad_samples = [s + sample_seq[3:20] for s in ["AAA", "ATA", "AGA", "ACG"]]

    regression_data['bad_sample_runs'] = []

    for bad_sample in bad_samples:
        result_bad = get_trigger_mfe_scores_by_risearch(
            bad_sample,
            name_to_seq,
            interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
            minimum_score=900,
            neighborhood=30,
            parsing_type='2',
            transpose=True
        )
        mfe_scores_bad = get_mfe_scores(result_bad, '2')

        # Append each iteration's result to the list
        regression_data['bad_sample_runs'].append({
            'input_sequence': bad_sample,
            'risearch_result': result_bad,
            'mfe_scores': mfe_scores_bad
        })

    # --- 5. Perform Regression Check ---
    # This checks the current `regression_data` dict against the stored YAML file.
    data_regression.check(regression_data)



def test_risearch_gfp_modified_original(data_regression):
    gfp_seq = get_gfp_first_exp(gap=0)
    sample_seq = gfp_seq[:20]
    name_to_seq = {f"gfp_seq": gfp_seq}

    # Initialize a dictionary to hold all data we want to snapshot
    regression_data = {}

    # Capture context data (replacing your initial print statements)
    regression_data['metadata'] = {
        'gfp_ontarget_slice': gfp_seq[695:714],
        'sample_seq': sample_seq,
        'sample_antisense': get_antisense(sample_seq)
    }

    # --- 3. Test "Good" Sample ---
    result_good = get_trigger_mfe_scores_by_risearch(
        sample_seq,
        name_to_seq,
        interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
        minimum_score=900,
        neighborhood=30,
        parsing_type='2',
        transpose=False
    )
    mfe_scores_good = get_mfe_scores(result_good, '2')

    print(mfe_scores_good)

    # Add results to our data snapshot
    regression_data['good_sample_run'] = {
        'risearch_result': result_good,
        'mfe_scores': mfe_scores_good
    }

    # --- 4. Test "Bad" Samples ---
    bad_samples = [s + sample_seq[3:20] for s in ["AAA", "ATA", "AGA", "ACG"]]

    regression_data['bad_sample_runs'] = []

    for bad_sample in bad_samples:
        result_bad = get_trigger_mfe_scores_by_risearch(
            bad_sample,
            name_to_seq,
            interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
            minimum_score=900,
            neighborhood=30,
            parsing_type='2',
            transpose=False
        )
        mfe_scores_bad = get_mfe_scores(result_bad, '2')

        # Append each iteration's result to the list
        regression_data['bad_sample_runs'].append({
            'input_sequence': bad_sample,
            'risearch_result': result_bad,
            'mfe_scores': mfe_scores_bad
        })

    # --- 5. Perform Regression Check ---
    # This checks the current `regression_data` dict against the stored YAML file.
    data_regression.check(regression_data)
