from .utils import *
from .data_handling import *
from ..features.mod_features import compute_mod_min_distance_to_3prime
from .populate.populate_sense_accessibility import populate_sense_accessibility
from .populate.populate_cai import populate_cai_for_aso_dataframe
from ..features.RNaseH_features import rnaseh1_dict, compute_rnaseh1_score

FLANK_SIZE = 120
ACCESS_SIZE = 13
SEED_SIZE = 13
SEED_SIZES = [SEED_SIZE * m for m in range(1, 4)]
ACCESS_WIN_SIZE = 80


def add_RNaseH1_Krel(df, exp='R7_krel'):
    best_window_start_krel = {
        'R4a_krel': {10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 3, 18: 2, 19: 4, 20: 3, 21: 0, 22: 0, 25: 0},
        'R4b_krel': {10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 1, 18: 3, 19: 1, 20: 3, 21: 0, 22: 0, 25: 0},
        'R7_krel': {10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 3, 17: 2, 18: 4, 19: 6, 20: 4, 21: 0, 22: 0, 25: 0},
    }
    weights = rnaseh1_dict(exp)

    def score_row(row):
        length = len(row['Sequence'])
        pos = best_window_start_krel.get(exp, {}).get(length, 0)
        return compute_rnaseh1_score(row['Sequence'], weights, window_start=pos)

    col_name = f"RNaseH1_Krel_score_{exp}"
    df[col_name] = df.apply(score_row, axis=1)
    return df


def fill_df(df, gene, genes_lst, gene_to_data, tp=24, vol=1000, mod_type="moe", mod_pattern=None):
    assert (len(genes_lst) == 1)
    if mod_pattern is None:
        mod_pattern = mod_type_dic[mod_type]
        df['mod_pattern'] = mod_pattern
    else:
        df['mod_pattern'] = ['d' * int(l) for l in df[SENSE_LENGTH]]

    df[CANONICAL_GENE] = gene
    df[CELL_LINE_ORGANISM] = 'human'
    df[INHIBITION] = 0  ## for function
    df = get_populated_df_with_structure_features(df, genes_lst, gene_to_data)
    df[TREATMENT_PERIOD] = tp  # keep constant for all
    df[VOLUME] = vol  # keep constant for all
    df['log_volume'] = np.log(df[VOLUME])
    df['normalized_start'] = df[SENSE_START] / len(gene_to_data[gene].full_mrna)
    df['normalized_sense_start_from_end'] = df['sense_start_from_end'] / len(gene_to_data[gene].full_mrna)
    easy_to_populate = ['at_skew', 'gc_content', 'gc_content_3_prime_5', 'gc_skew', 'hairpin_score',
                        'homooligo_count', 'internal_fold', 'nucleotide_diversity', 'self_energy', 'stop_codon_count',
                        'at_rich_region_score', 'poly_pyrimidine_stretch']
    populate_features(df, easy_to_populate)
    fold_variants = [(40, 15)]
    df = get_populate_fold(df, genes_lst, gene_to_data, fold_variants=fold_variants)
    # TO DO:
    df.loc[:, 'Modification_min_distance_to_3prime'] = compute_mod_min_distance_to_3prime(
        mod_pattern)  ## gen more mod_patterns

    # df.loc[:, 'sense_avg_accessibility'] = 0.8526221963673779 # TODO - replace with the calculation
    populate_sense_accessibility(df, gene_to_data[gene])
    populate_cai_for_aso_dataframe(df, gene_to_data[gene])
    # df.loc[:, SENSE_AVG_ACCESSIBILITY] = 0.8526221963673779 # TODO - replace with the calculation
    # df.loc[:, 'CAI_score_global_CDS'] = 0.8526221963673779 # TODO - replace with the calculation
    # df.loc[:,'sense_avg_accessibility'] = 0.5 # takes long time to calc

    add_RNaseH1_Krel(df)
    return df


def fill_dfs(genes_lst, gene_to_data, mod_type='moe', tp=24, vol=1000):
    dfs = {}
    for gene in genes_lst:
        gene_info = gene_to_data[gene]
        dfs[gene] = get_init_df(gene_info.full_mrna, gene_info.exon_indices[-1][1] - gene_info.cds_start,
                                mod_type=mod_type)
    for gene, df in dfs.items():
        dfs[gene] = fill_df(df, gene, genes_lst, gene_to_data, mod_type=mod_type)

    return dfs

# def compute_sense_accessibility(row, flank_size, access_win_size, seed_sizes, access_size, min_gc=0, max_gc=100, gc_ranges=1):


#     try:
#         # Skip invalid rows
#         if row['sense_start'] == -1 or pd.isna(row['sense_with_flank_120nt']) or row['sense_with_flank_120nt'] == "":
#             return None

#         seq = row[f'sense_with_flank_{flank_size}nt']
#         sense_start = row['sense_start']
#         sense_length = row['sense_length']

#         # Calculate accessibility
#         df_access = AccessCalculator.calc(
#             seq, access_size,
#             min_gc, max_gc, gc_ranges,
#             access_win_size, seed_sizes
#         )

#         flank_start = max(0, sense_start - flank_size)
#         sense_start_in_flank = sense_start - flank_start
#         sense_end_in_flank = sense_start_in_flank + sense_length

#         if 0 <= sense_start_in_flank < len(df_access) and sense_end_in_flank <= len(df_access):
#             values = df_access['avg_access'].iloc[sense_start_in_flank:sense_end_in_flank].dropna()
#             return values.mean() if not values.empty else None
#         else:
#             return None

#     except Exception as e:
#         print(f"Error at row {row.name} | seq start: {row['sense_start']} | error: {e}")
#         return None


# for gene, df in dfs.items():
#     FLANKED_SENSE_COL = f'sense_with_flank_{FLANK_SIZE}nt'

#     val = gene_to_data[gene].full_mrna
#     df['pre_mrna_sequence'] = [val] * len(df)


#     # Create new column with flanked sequences
#     df[FLANKED_SENSE_COL] = df.apply(
#     lambda row: get_sense_with_flanks(
#         row['pre_mrna_sequence'],
#         row['sense_start'],
#         row['sense_length'],
#         flank_size=FLANK_SIZE
#     ) if row['sense_start'] != -1 else "",  # Handle cases where sense was not found
#     axis=1
#     )


#     batch_size = 500
#     for start_idx in range(0, len(df), batch_size):
#         end_idx = min(start_idx + batch_size, len(df))
#         batch = df.iloc[start_idx:end_idx].copy()

#         print(f"Processing rows {start_idx} to {end_idx}...")

#         batch['sense_avg_accessibility'] = batch.apply(
#             compute_sense_accessibility,
#             axis=1,
#             flank_size=FLANK_SIZE,
#             access_win_size=ACCESS_WIN_SIZE,
#             seed_sizes=SEED_SIZES,
#             access_size=ACCESS_SIZE,
#         )

#         # Save batch to the new folder
#         batch.to_csv(f"out/batch_{start_idx}_{end_idx}.csv", index=False)
