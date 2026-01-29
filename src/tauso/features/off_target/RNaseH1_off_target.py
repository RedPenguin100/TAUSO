import pandas as pd
import numpy as np
from .off_target_functions import parse_risearch_output
from ...hybridization.fast_hybridization import get_trigger_mfe_scores_by_risearch, Interaction
from ...new_model.consts_dataframe import SEQUENCE
from ...util import get_antisense


def off_target_specific_seq(aso_df, gene_name, gene_to_data, cutoff):
    sequence = gene_to_data[gene_name].full_mrna

    name_to_seq = {gene_name: sequence}
    RT = 0.616
    results = []

    feature_name = f"off_target_single_{gene_name}_c{cutoff}"

    for idx_val, aso_seq in zip(aso_df['index'], aso_df[SEQUENCE]):
        trigger = get_antisense(aso_seq)

        result_dict = get_trigger_mfe_scores_by_risearch(
            trigger,
            name_to_seq,
            minimum_score=cutoff,
            parsing_type='2',
            interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
            transpose=True
        )

        result_df = parse_risearch_output(result_dict)

        # Handle cases where result_df might be empty (no hits above cutoff)
        if not result_df.empty and 'energy' in result_df.columns:
            off_target_score = np.exp(-RT * result_df['energy']).sum()
        else:
            off_target_score = 0.0

        results.append({"index": idx_val, feature_name: off_target_score})

    score_df = pd.DataFrame(results)

    return score_df, feature_name
