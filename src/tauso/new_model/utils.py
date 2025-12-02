# Check the wellness of fit
from xgboost import XGBRanker
import pickle
from ..util import get_antisense
import pandas as pd
from ..genome.read_human_genome import get_locus_to_data_dict
from pathlib import Path

mod_type_dic = {'moe': 'MMMMMddddddddddMMMMM', 'lna': 'LLLddddddddddLLL'}

metric = 'correct_log_inhibition2'

SEQUENCE = 'Sequence'
INHIBITION = 'Inhibition(%)'
CANONICAL_GENE = 'Canonical Gene Name'
CELL_LINE_ORGANISM = 'Cell line organism'
VOLUME = 'ASO_volume(nM)'
CHEMICAL_PATTERN = 'Chemical_Pattern'
TREATMENT_PERIOD = 'Treatment_Period(hours)'
CELL_LINE = 'Cell_line'
TRANSFECTION = 'Transfection'
DENSITY = 'Density(cells/well)'
DENSITY_UPDATED = 'Density(cells_per_well)'  # Avoiding /
MODIFICATION = 'Modification'
PREMRNA_FOUND = 'pre_mrna_found'
SENSE_START = 'sense_start'
SENSE_LENGTH = 'sense_length'
SENSE_TYPE = 'sense_type'

parent = Path(__file__).parent

selected_features = [
    TREATMENT_PERIOD,
    'at_skew',
    'CAI_score_global_CDS',
    'stop_codon_count',
    'sense_avg_accessibility',
    'on_target_fold_openness_normalized40_15',
    'sense_utr',
    'nucleotide_diversity',
    'internal_fold',
    'normalized_start',
    'RNaseH1_Krel_score_R7_krel',  # renamed to best
    'hairpin_score',
    'Modification_min_distance_to_3prime',
    'at_rich_region_score'
]

tbl = str.maketrans("ACGTUacgtuNn", "TGCAAtgcaaNn")


def get_init_df(target_mrna, end, mod_type='moe'):
    s = len(mod_type_dic[mod_type])
    candidates = []
    sense_starts = []
    sense_lengths = []
    sense_starts_from_end = []
    set_candidates = set()
    counter = 0
    for i in range(0, len(target_mrna) - (s - 1)):
        target = target_mrna[i: i + s]
        if target in set_candidates:
            counter = + 1
            continue
        set_candidates.add(target)
        candidates.append(get_antisense(str(target)).upper().replace('T', 'U'))
        sense_starts.append(i)
        sense_lengths.append(s)
        sense_starts_from_end.append(end - i)
    df = pd.DataFrame(
        {SEQUENCE: candidates, SENSE_START: sense_starts,
         SENSE_LENGTH: sense_lengths, "sense_start_from_end": sense_starts_from_end})
    return df


def create_gene_to_data(genes_u):  # can work only for our known list
    cache_path = Path("./gene_to_data_simple_cache.pickle")
    if not cache_path.exists():
        gene_to_data = get_locus_to_data_dict(include_introns=True, gene_subset=genes_u)

    else:
        print("pickle exist")
        with open(cache_path, 'rb') as f:
            gene_to_data = pickle.load(f)
    return gene_to_data


def create_and_load_model(json_weight=str(parent / "model.json"), seed=42):
    model = XGBRanker(objective='rank:ndcg', ndcg_exp_gain=False, lambdarank_pair_method="topk",
                      lambdarank_num_pair_per_sample=200,
                      seed=seed, n_jobs=-1)
    model.load_model(json_weight)
    return model

# def plot_res(model,X,y,df_copy):
#     predicted = model.predict(X[selected_features].values)

#     corr, _ = pearsonr(predicted, y)
#     corrs, _ = spearmanr(predicted, y)

#     print(
#         f"Pearson {corr}, Spearman {corrs}")


#     y_true_corrected = df_copy[metric + 'ndcg'].to_numpy()

#     # ---- New metrics ----
#     # NDCG@50
#     ndcg50 = ndcg_score(y_true_corrected.reshape(1, -1), predicted.reshape(1, -1), k=200)
#     # NDCG (all, just drop k)
#     ndcg_all = ndcg_score(y_true_corrected.reshape(1, -1), predicted.reshape(1, -1))
#     # Precision@50 (manual overlap of top-50 by truth vs pred)
#     precisions = []
#     values = [50, 100]
#     for K in values:
#         pred_top_idx = np.argpartition(predicted, -K)[-K:]
#         true_top_idx = np.argpartition(y_true_corrected, -K)[-K:]
#         precisions.append(len(set(pred_top_idx) & set(true_top_idx)) / K)

#     print(f"NDCG@50: {ndcg50:.4f}, NDCG(all): {ndcg_all:.4f}, Precision: {values, precisions}")
