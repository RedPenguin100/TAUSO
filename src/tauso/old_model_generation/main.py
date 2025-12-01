from pathlib import Path
from typing import Optional, List, Dict, Any, Union, Tuple

import numpy as np
import pandas as pd

from xgboost import XGBRanker

from tauso.features.mod_features import compute_mod_min_distance_to_3prime
from tauso.genome.read_human_genome import get_locus_to_data_dict
from tauso.new_model.data_handling import get_populate_fold, populate_features, get_populated_df_with_structure_features
from tauso.new_model.feature_creation import add_RNaseH1_Krel
from tauso.new_model.populate.populate_cai import populate_cai_for_aso_dataframe
from tauso.new_model.populate.populate_sense_accessibility import populate_sense_accessibility
from tauso.off_target.search import find_all_gene_off_targets
from tauso.old_model_generation.consts_dataframe import SEQUENCE, SENSE_START, SENSE_LENGTH, MOD_TYPE_DICT, \
    TREATMENT_PERIOD, VOLUME, INHIBITION, CELL_LINE_ORGANISM, CANONICAL_GENE
from tauso.util import get_antisense


def get_target_sequence(
        target_name: str,
        custom_sequence: Optional[str] = None
) -> str:
    """
    Retrieves the mRNA sequence either from the custom input or
    by querying the local genome database.
    """
    if custom_sequence:
        return custom_sequence.upper().replace('U', 'T')

    # Use the optimized subset loader from read_human_genome.py
    print(f"Fetching sequence for {target_name}...")
    locus_data = get_locus_to_data_dict(include_introns=False, gene_subset=[target_name])

    if target_name not in locus_data:
        raise ValueError(f"Gene {target_name} not found in genome database.")

    # Extract the full mRNA sequence from the LocusInfo object
    return locus_data[target_name].full_mrna


DEFAULT_ASO_LENGTH = 20


def get_init_df(target_mrna: str, mod_type: str = 'moe') -> pd.DataFrame:
    """
    Generates ASO candidates by tiling the mRNA.

    Changes from legacy:
    - Removed 'end' parameter.
    - uses len(target_mrna) to calculate 'sense_start_from_end'.
    """
    # Determine ASO length (s) based on chemistry
    if mod_type in MOD_TYPE_DICT:
        s = len(MOD_TYPE_DICT[mod_type])
    else:
        s = 20  # Default fallback

    target_len = len(target_mrna)

    candidates = []
    sense_starts = []
    sense_lengths = []
    sense_starts_from_end = []

    # Deduplication set
    set_candidates = set()

    # Tile the sequence
    for i in range(0, target_len - (s - 1)):
        target = target_mrna[i: i + s]

        # Deduplicate based on the SENSE target sequence
        if target in set_candidates:
            continue
        set_candidates.add(target)

        # Store Data
        candidates.append(get_antisense(str(target)))
        sense_starts.append(i)
        sense_lengths.append(s)

        # distance from 3' end (assuming full_mrna is the context)
        # Note: Index 0 is 5' end. Distance from 3' end decreases as i increases.
        sense_starts_from_end.append(target_len - i)

    df = pd.DataFrame({
        SEQUENCE: candidates,
        SENSE_START: sense_starts,
        SENSE_LENGTH: sense_lengths,
        "sense_start_from_end": sense_starts_from_end
    })

    return df


def enrich_shared_features(
        df: pd.DataFrame,
        gene: str,
        genes_lst: List[str],
        gene_to_data: Dict[str, Any],
        SEQUENCES: Dict[str, str],
        tp: float,
        vol: float
) -> pd.DataFrame:
    """
    CALCULATED ONCE:
    Computes expensive features that depend only on the sequence or gene,
    not the specific modification pattern string.
    """
    # 1. Constants
    df[CANONICAL_GENE] = gene
    df[CELL_LINE_ORGANISM] = 'human'
    df[TREATMENT_PERIOD] = tp
    df[VOLUME] = vol
    df['log_volume'] = np.log(df[VOLUME])
    df[INHIBITION] = 0

    # 2. Normalization
    seq_len = len(SEQUENCES[gene])
    df['normalized_start'] = df[SENSE_START] / seq_len
    df['normalized_sense_start_from_end'] = df['sense_start_from_end'] / seq_len

    # 3. Structural Features (EXPENSIVE)
    # Folding the candidate and checking target structure
    df = get_populated_df_with_structure_features(df, genes_lst, gene_to_data)

    # 4. Sequence Features (GC, Skew, etc.)
    easy_to_populate = [
        'at_skew', 'gc_content', 'gc_content_3_prime_5', 'gc_skew', 'hairpin_score',
        'homooligo_count', 'internal_fold', 'nucleotide_diversity', 'self_energy',
        'stop_codon_count', 'at_rich_region_score', 'poly_pyrimidine_stretch'
    ]
    populate_features(df, easy_to_populate)

    # 5. Fold Variants
    fold_variants = [(40, 15)]
    df = get_populate_fold(df, genes_lst, gene_to_data, fold_variants=fold_variants)

    # 6. Gene-Level Features (EXPENSIVE)
    # Accessibility and CAI depend on the target gene, not the ASO chemistry
    populate_sense_accessibility(df, gene_to_data[gene])
    populate_cai_for_aso_dataframe(df, gene_to_data[gene])

    return df


def enrich_chemistry_features(df: pd.DataFrame, mod_type: str) -> pd.DataFrame:
    """
    CALCULATED PER CHEMISTRY:
    Adds lightweight tags and features specific to the modification pattern.
    """
    mod_pattern = MOD_TYPE_DICT[mod_type]
    df['mod_pattern'] = mod_pattern

    # Distance depends on the specific pattern string
    df['Modification_min_distance_to_3prime'] = compute_mod_min_distance_to_3prime(mod_pattern)

    # RNase H1 score might depend on the gapmer pattern
    add_RNaseH1_Krel(df)

    return df


parent = Path(__file__).parent


def create_and_load_model(json_weight=str(parent / "model.json"), seed=42):
    model = XGBRanker(objective='rank:ndcg', ndcg_exp_gain=False, lambdarank_pair_method="topk",
                      lambdarank_num_pair_per_sample=200,
                      seed=seed, n_jobs=-1)
    model.load_model(json_weight)
    return model


def enrich_with_off_targets(
        df: pd.DataFrame,
        genome: str = 'GRCh38',
        max_mismatches: int = 3
) -> pd.DataFrame:
    """
    Runs the new Bowtie-based off-target search (tauso.off_target.search)
    Only runs on the provided rows (efficient for top-k).
    """
    print(f"Running off-target analysis for {len(df)} candidates...")

    results = []

    for _, row in df.iterrows():
        aso_seq = row[SEQUENCE]

        # Call the new search.py function
        try:
            hits_df = find_all_gene_off_targets(
                aso_seq,
                genome=genome,
                max_mismatches=max_mismatches
            )

            if hits_df.empty:
                count_0mm = 0
                count_total = 0
                summary = "No off-targets"
            else:
                count_0mm = len(hits_df[hits_df['mismatches'] == 0])
                count_total = len(hits_df)

                # Create a concise summary string (e.g., "GAPDH(0mm); ACTB(1mm)")
                # Prioritize low mismatches and genes
                top_hits = hits_df.sort_values(['mismatches', 'gene_name']).head(3)
                summary_parts = []
                for _, h in top_hits.iterrows():
                    gname = h['gene_name'] if h['gene_name'] else h['chrom']
                    summary_parts.append(f"{gname}({h['mismatches']})")
                summary = "; ".join(summary_parts)

        except Exception as e:
            print(f"Warning: Search failed for {aso_seq}: {e}")
            count_0mm = -1
            count_total = -1
            summary = "Search Error"

        results.append({
            'off_target_0mm': count_0mm,
            'off_target_total': count_total,
            'off_target_summary': summary
        })

    ot_df = pd.DataFrame(results)
    # Join back to original DF
    return pd.concat([df.reset_index(drop=True), ot_df], axis=1)


DEFAULT_GENOME = 'GRCh38'


def design_asos(
        gene_name: str,
        organism: str = 'human',
        custom_sequence: Optional[str] = None,
        top_k: int = 10,
        include_details: bool = True,
        run_off_target: bool = True,
        return_seq: bool = False,
        mod_type: List[str] = ['moe', 'lna'],

) -> Union[pd.DataFrame, Tuple[pd.DataFrame, str]]:
    """
    Main pipeline entry point.

    Args:
        gene_name: Name of gene (e.g. 'MALAT1')
        organism: 'human' triggers DB lookup and off-target search.
        custom_sequence: Optional override for sequence.
        top_k: Number of candidates to return.
        include_details: If True, returns feature columns.
        run_off_target: If True, runs Bowtie search on top candidates.
        mod_type: List of chemistries to score ('moe', 'lna').

    Returns:
        DataFrame of designed ASOs.
    """

    # 1. Get Sequence
    target_seq = get_target_sequence(gene_name, custom_sequence)

    if not gene_name:
        raise ValueError("TODO")
    locus_data = get_locus_to_data_dict(include_introns=False, gene_subset=[gene_name])

    # 1. Group chemistries by ASO length
    # Example: {20: ['moe', 'gapmer'], 16: ['lna']}
    chem_by_len = {}
    for chem in mod_type:
        # Determine length based on pattern
        pattern_len = len(MOD_TYPE_DICT.get(chem, 'A' * 20))
        if pattern_len not in chem_by_len:
            chem_by_len[pattern_len] = []
        chem_by_len[pattern_len].append(chem)

    results_list = []

    # 2. Iterate by Length Group
    for length, chem_group in chem_by_len.items():
        print(f"Processing length {length} for chemistries: {chem_group}")

        # A. Generate Candidates ONCE for this length
        # Note: We pass a dummy mod_type just to trigger the correct length in get_init_df
        # or we explicitly handle length in get_init_df if you updated it to take 'length' arg.
        # Using the first chem in group to set length context:
        base_df = get_init_df(target_seq, mod_type=chem_group[0])

        # B. Calculate Expensive Features ONCE for this DataFrame
        print(f"  - Pre-calculating heavy features for {len(base_df)} candidates...")
        try:
            # We pass a dictionary for SEQUENCES as expected by legacy functions
            rich_df = enrich_shared_features(
                base_df,
                gene=gene_name,
                genes_lst=[gene_name],
                gene_to_data=locus_data,
                SEQUENCES={gene_name: target_seq},
                tp=24,
                vol=1000
            )
        except Exception as e:
            print(f"Feature calculation failed: {e}")
            continue

        # C. Apply specific chemistry tags and Predict
        for chem in chem_group:
            print(f"  - Scoring {chem}...")
            try:
                # Copy the rich dataframe (cheap compared to recalculating features)
                chem_df = rich_df.copy()

                # Add lightweight chemistry columns
                chem_df = enrich_chemistry_features(chem_df, mod_type=chem)

                # Predict
                # (We inline the prediction part of run_scoring_model here)
                selected_features = [
                    TREATMENT_PERIOD, 'at_skew', 'CAI_score_global_CDS', 'stop_codon_count',
                    'sense_avg_accessibility', 'on_target_fold_openness_normalized40_15',
                    'sense_utr', 'nucleotide_diversity', 'internal_fold', 'normalized_start',
                    'RNaseH1_Krel_score_R7_krel', 'hairpin_score',
                    'Modification_min_distance_to_3prime', 'at_rich_region_score'
                ]

                # Ensure columns exist
                for col in selected_features:
                    if col not in chem_df.columns: chem_df[col] = 0.0

                model = create_and_load_model()
                chem_df['score'] = model.predict(chem_df[selected_features].values)

                results_list.append(chem_df)

            except Exception as e:
                print(f"Error scoring {chem}: {e}")
                continue

    if not results_list:
        raise RuntimeError("No candidates could be scored.")

    full_results = pd.concat(results_list)

    # 4. Sort and Filter Top K
    # Sort by score descending
    full_results = full_results.sort_values(by='score', ascending=False)

    # Take top K *before* expensive off-target search
    top_results = full_results.head(top_k).copy()

    # 5. Run Off-Target Analysis (Expensive Step)
    if run_off_target and organism == 'human':
        top_results = enrich_with_off_targets(top_results, genome=DEFAULT_GENOME)

    # 6. Select Columns for Output
    # Basic identification info
    output_cols = [SEQUENCE, SENSE_START, 'mod_pattern', 'score', 'gc_content']

    # Off-target info
    if run_off_target and organism == 'human':
        output_cols.extend(['off_target_0mm', 'off_target_summary'])

    # Detailed features if requested
    if include_details:
        # Add all other columns present in df
        all_cols = list(top_results.columns)
        # Ensure we don't duplicate
        extra_cols = [c for c in all_cols if c not in output_cols]
        output_cols.extend(extra_cols)

    if return_seq:
        return top_results[output_cols], target_seq

    return top_results[output_cols]


# --- Demo Execution ---
if __name__ == "__main__":
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 1000)

    print("Running Design Pipeline Demo...")
    try:
        # This requires the genome DB to be set up via cli.py
        results = design_asos(
            gene_name="DDX11L1",
            top_k=5,
            run_off_target=True  # This triggers the new Bowtie search
        )
        print(results.columns)

        print("\nTop ASO Candidates:")
        print(results[['Sequence', 'mod_pattern', 'score', 'sense_avg_accessibility', 'off_target_total']])
    except Exception as e:
        print(f"\nDemo failed (likely due to missing DB or Models): {e}")
