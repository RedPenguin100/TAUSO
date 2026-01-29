from ...hybridization.fast_hybridization import get_trigger_mfe_scores_by_risearch, Interaction
from .off_target_functions import parse_risearch_output, aggregate_off_targets

from ...new_model.consts_dataframe import SEQUENCE, CANONICAL_GENE
from ...util import get_antisense

import numpy as np



class AggregationMethod:
    ARTM_log = "ARTM_log"
    ARTM = "ARTM"
    ARTM_weighted = "ARTM_weighted"
    GEO = "GEO"
    RANKED = "RANKED"
    MECH = "MECH"


def compute_general_vector_and_cache(ASO_df, gene_to_data, general_df, cutoff, method):
    """
    Calculates the General Vector for all rows and caches the RIsearch energies.
    """
    index_scores = {}
    energy_cache = {}  # Key: ASO_Seq, Value: { Gene: (Energy, Sequence) }

    # 1. Prepare General Data Maps
    general_seq_map = {}
    general_exp_map = {}

    norm_col = next((c for c in general_df.columns if 'expression_norm' in c), 'expression_norm')
    for _, row in general_df.iterrows():
        gene = row['Gene'].split()[0]
        seq = gene_to_data[gene].full_mrna

        general_seq_map[gene] = seq
        general_exp_map[gene] = (row.get('expression_TPM', 0), row.get(norm_col, row.get('expression_norm', 0)))

    # 2. Iterate Unique ASO Sequences (Avoid re-calculating identical ASOs)
    unique_pairs = ASO_df[[SEQUENCE, CANONICAL_GENE]].drop_duplicates()

    for _, row in unique_pairs.iterrows():
        aso_seq = row[SEQUENCE]
        target_gene = row[CANONICAL_GENE]  # Capture the specific gene for this pair

        trigger = get_antisense(aso_seq)

        # --- RIsearch (Heavy Lifting) ---
        result_dict = get_trigger_mfe_scores_by_risearch(
            trigger,
            general_seq_map,
            minimum_score=cutoff,
            interaction_type=Interaction.RNA_DNA_NO_WOBBLE,
            parsing_type='2',
            transpose=True
        )
        result_df = parse_risearch_output(result_dict)
        result_df_agg = aggregate_off_targets(result_df)

        # Cache Energies AND Sequences
        # Key is now a tuple (Sequence, Gene) to maintain uniqueness in the cache
        gene_energy_seq = {}
        if not result_df_agg.empty:
            for _, r in result_df_agg.iterrows():
                g = r['target']
                e = r['energy']
                gene_energy_seq[g] = e

        energy_cache[(aso_seq, target_gene)] = gene_energy_seq

        # Calculate General Score
        simple_energies = {g: e for g, e in gene_energy_seq.items()}
        score = calculate_score_helper(simple_energies, general_exp_map, method)

        # Assign this score to rows matching BOTH Sequence AND Canonical Gene
        # Note: We use the variables SEQUENCE and CANONICAL_GENE assuming they hold the column names
        mask = (ASO_df[SEQUENCE] == aso_seq) & (ASO_df[CANONICAL_GENE] == target_gene)
        relevant_indices = ASO_df[mask]['index']

        for idx in relevant_indices:
            index_scores[idx] = score

    return index_scores, energy_cache


def compute_specific_vector_optimized(specific_df, gene_to_data, ASO_df, energy_cache, cutoff, method):
    """
    Calculates scores for a specific cell line.
    OPTIMIZATION: Reuses cached energies if the gene sequence is identical to General.
    """
    index_scores = {}

    # 1. Prepare Specific Data
    specific_exp_map = {}
    specific_seq_map = {}

    norm_col = next((c for c in specific_df.columns if 'expression_norm' in c), 'expression_norm')

    for _, row in specific_df.iterrows():
        gene = row['Gene'].split()[0]
        seq = gene_to_data[gene].full_mrna

        specific_seq_map[gene] = seq
        specific_exp_map[gene] = (row.get('expression_TPM', 0), row.get(norm_col, row.get('expression_norm', 0)))

    # 2. Iterate ASOs
    for _, row in ASO_df.iterrows():
        index = row['index']
        aso_seq = row[SEQUENCE]
        target_gene = row[CANONICAL_GENE]

        # Retrieve Cache: {Gene: (Energy, General_Seq)}
        cached_data = energy_cache.get(aso_seq, {})

        # 3. Identify "Delta" Genes (Mutations or New Genes)
        delta_seq_map = {}
        final_energies = {}

        for gene, spec_seq in specific_seq_map.items():
            if gene == target_gene: continue

            # Check Cache
            if gene in cached_data:
                gen_energy = cached_data[gene]
                # MATCH: Reuse cached energy
                final_energies[gene] = gen_energy
                continue

            # MISMATCH: Sequence differs (or new gene) -> Add to Delta
            delta_seq_map[gene] = spec_seq

        # 4. Run RIsearch ONLY on Delta Genes
        if delta_seq_map:
            trigger = get_antisense(aso_seq)
            delta_result = get_trigger_mfe_scores_by_risearch(
                trigger,
                delta_seq_map,
                minimum_score=cutoff,
                parsing_type='2'
            )
            delta_df = parse_risearch_output(delta_result)
            delta_agg = aggregate_off_targets(delta_df)

            # Add delta energies to final map
            if not delta_agg.empty:
                for _, r in delta_agg.iterrows():
                    final_energies[r['target']] = r['energy']

        # 5. Calculate Final Score
        score = calculate_score_helper(final_energies, specific_exp_map, method)
        index_scores[index] = score

    return index_scores


def calculate_score_helper(energy_dict, expression_dict, method):
    """
    Standardizes the scoring math to avoid code duplication.
    """
    score = 0.0
    RT = 0.616

    if not energy_dict:
        return 0.0

    # Pre-filter for valid interactions (Energy < 0) and Expression > 0
    # This speeds up the loops significantly
    valid_targets = []
    for gene, energy in energy_dict.items():
        exp = expression_dict.get(gene, 0)[0]
        # Use 1e6 scaling for TPM-based methods (1, 2, 3, 4, 5)
        # Method 0 uses norm, which shouldn't be scaled by 1e6
        if method != 0:
            exp = exp / 1e6

        if energy < 0 and exp > 0:
            valid_targets.append((gene, energy, exp))

    if not valid_targets:
        return 0.0

    if method == AggregationMethod.ARTM_log:  # Normalized Arithmetic (Direct exp, no 1e6 scaling)
        # Re-loop because we scaled exp above, but method 0 usually takes raw norm
        score = 0.0
        for gene, energy in energy_dict.items():
            exp = expression_dict.get(gene, 0)[1]  # Raw Norm
            if energy < 0 and exp > 0:
                score += energy * exp

    elif method == AggregationMethod.ARTM:  # TPM Arithmetic
        for _, energy, exp in valid_targets:
            score += energy * exp

    elif method == AggregationMethod.ARTM_weighted:  # Weighted Arithmetic (e^2)
        for _, energy, exp in valid_targets:
            score += energy * (exp ** 2)

    elif method == AggregationMethod.GEO:  # Geometric
        log_sum = 0.0
        for _, energy, exp in valid_targets:
            log_sum += exp * np.log(-energy)
        score = log_sum

    elif method == AggregationMethod.RANKED:  # Ranked
        # Sort by expression descending
        valid_targets.sort(key=lambda x: x[2], reverse=True)

        for rank, (_, energy, exp) in enumerate(valid_targets, start=1):
            batch = (rank - 1) // 10
            weight = 10 / (2 ** batch)
            score += energy * exp * weight

    elif method == AggregationMethod.MECH:  # Mechanical Statistics
        for _, energy, exp in valid_targets:
            score += exp * np.exp(-energy / RT)


    # Normalize Mechanical Statistics if required (usually done outside, but can be done here if context allows)
    # Note: Normalization requires knowing the max of the whole vector, so strictly speaking
    # method 5 normalization should happen after collecting all scores.
    # For now, this returns the raw Sum.

    return score