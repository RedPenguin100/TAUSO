import primer3
import numpy as np

from ViennaRNA import RNA

from ..features.vienna_fold import calculate_energies, get_weighted_energy
from ..util import get_antisense
from ..features.seq_features import palindromic_fraction, homooligo_count, hairpin_score, seq_entropy, \
    gc_skew, at_skew, \
    nucleotide_diversity, stop_codon_count, get_gc_content, at_rich_region_score, poly_pyrimidine_stretch
from .consts_dataframe import CELL_LINE_ORGANISM, CANONICAL_GENE, SEQUENCE, SENSE_LENGTH, SENSE_START
from .utils import INHIBITION
from ..genome.read_human_genome import get_locus_to_data_dict


def get_unique_human_genes(all_data):
    all_data_human = all_data[all_data[CELL_LINE_ORGANISM] == 'human']
    all_data_human_no_nan = all_data_human.dropna(subset=[INHIBITION]).copy()

    genes = all_data_human_no_nan[CANONICAL_GENE].copy()
    genes_u = list(set(genes))

    genes_u.remove('HBV')
    genes_u.remove('negative_control')

    return genes_u


def get_gene_to_data(genes_u):
    import pickle
    from ..consts import CACHE_DIR

    cache_path = CACHE_DIR / 'gene_to_data_simple_cache.pickle'

    # TODO: hash the pickled file to avoid mis-reads
    if not cache_path.exists():
        gene_to_data = get_locus_to_data_dict(include_introns=True, gene_subset=genes_u)
        with open(cache_path, 'wb') as f:
            pickle.dump(gene_to_data, f)
    else:
        with open(cache_path, 'rb') as f:
            gene_to_data = pickle.load(f)

    return gene_to_data


def get_populated_df_with_structure_features(df, genes_u, gene_to_data):
    """
    Populate "the data" df with features like exon/intron, start of the sense strand, if found.
    """
    df_copy = df.copy()
    all_data_human = df_copy[df_copy[CELL_LINE_ORGANISM] == 'human']
    all_data_human_no_nan = all_data_human.dropna(subset=[INHIBITION]).copy()
    all_data_human_gene = all_data_human_no_nan[all_data_human_no_nan[CANONICAL_GENE].isin(genes_u)].copy()
    SENSE_START = 'sense_start'
    SENSE_START_FROM_END = 'sense_start_from_end'
    SENSE_LENGTH = 'sense_length'
    SENSE_TYPE = 'sense_type'
    SENSE_EXON = 'sense_exon'
    SENSE_INTRON = 'sense_intron'
    SENSE_UTR = 'sense_utr'

    found = 0
    all_data_human_gene[SENSE_START] = np.zeros_like(all_data_human_gene[CANONICAL_GENE], dtype=int)
    all_data_human_gene[SENSE_START_FROM_END] = np.zeros_like(all_data_human_gene[CANONICAL_GENE], dtype=int)
    all_data_human_gene[SENSE_LENGTH] = np.zeros_like(all_data_human_gene[CANONICAL_GENE], dtype=int)
    all_data_human_gene[SENSE_EXON] = np.zeros_like(all_data_human_gene[CANONICAL_GENE], dtype=int)
    all_data_human_gene[SENSE_INTRON] = np.zeros_like(all_data_human_gene[CANONICAL_GENE], dtype=int)
    all_data_human_gene[SENSE_UTR] = np.zeros_like(all_data_human_gene[CANONICAL_GENE], dtype=int)
    all_data_human_gene[SENSE_TYPE] = "NA"
    for index, row in all_data_human_gene.iterrows():
        gene_name = row[CANONICAL_GENE]
        locus_info = gene_to_data[gene_name]
        pre_mrna = locus_info.full_mrna
        antisense = row[SEQUENCE]
        sense = get_antisense(antisense)
        idx = pre_mrna.find(sense)
        all_data_human_gene.loc[index, SENSE_START] = idx
        all_data_human_gene.loc[index, SENSE_START_FROM_END] = np.abs(
            locus_info.exon_indices[-1][1] - locus_info.cds_start - idx)
        all_data_human_gene.loc[index, SENSE_LENGTH] = len(antisense)
        if idx != -1:
            genome_corrected_index = idx + locus_info.cds_start
            found = False
            for exon_indices in locus_info.exon_indices:
                # print(exon[0], exon[1])
                if exon_indices[0] <= genome_corrected_index <= exon_indices[1]:
                    all_data_human_gene.loc[index, SENSE_TYPE] = 'exon'
                    all_data_human_gene.loc[index, SENSE_EXON] = 1
                    found = True
                    break
            for intron_indices in locus_info.intron_indices:
                # print(exon[0], exon[1])
                if intron_indices[0] <= genome_corrected_index <= intron_indices[1]:
                    all_data_human_gene.loc[index, SENSE_TYPE] = 'intron'
                    all_data_human_gene.loc[index, SENSE_INTRON] = 1
                    found = True
                    break
            for i, utr_indices in enumerate(locus_info.utr_indices):
                if utr_indices[0] <= genome_corrected_index <= utr_indices[1]:
                    all_data_human_gene.loc[index, SENSE_TYPE] = 'utr'
                    all_data_human_gene.loc[index, SENSE_UTR] = 1

                    found = True
                    break
        if not found:
            all_data_human_gene.loc[index, SENSE_TYPE] = 'intron'
    return all_data_human_gene


def get_populate_fold(df, genes_u, gene_to_data, fold_variants=[(40, 15)]):
    all_data_human_gene_premrna_no_nan = df.copy()

    # Comment out the long cases for quick running
    for (window_size, step_size) in fold_variants:

        on_target_fold = 'on_target_fold_openness' + str(window_size) + '_' + str(step_size)
        on_target_fold_normalized = 'on_target_fold_openness_normalized' + str(window_size) + '_' + str(step_size)
        all_data_human_gene_premrna_no_nan[on_target_fold] = np.zeros_like(all_data_human_gene_premrna_no_nan[SEQUENCE],
                                                                           dtype=np.float64)
        all_data_human_gene_premrna_no_nan[on_target_fold_normalized] = np.zeros_like(
            all_data_human_gene_premrna_no_nan[SEQUENCE], dtype=np.float64)

        for gene in genes_u:

            target = gene_to_data[gene].full_mrna
            gene_rows = all_data_human_gene_premrna_no_nan[all_data_human_gene_premrna_no_nan[CANONICAL_GENE] == gene]
            energies = calculate_energies(str(target), step_size, window_size)

            for index, row in gene_rows.iterrows():
                antisense = row[SEQUENCE]
                sense = get_antisense(antisense)
                l = row[SENSE_LENGTH]
                sense_start = row[SENSE_START]
                print("Sense start internal: ", sense_start)
                mean_fold = get_weighted_energy(sense_start, l, step_size, energies, window_size)
                mean_fold_end = get_weighted_energy(sense_start, l, step_size, energies, window_size)
                mean_fold_start = get_weighted_energy(sense_start, l, step_size, energies, window_size)
                if mean_fold > 100:
                    print(energies)
                    print("Weird: ", mean_fold)
                    print("Sense_start ", sense_start)
                    print("Sense_length ", l)
                    print("Gene: ", gene)
                    mean_fold = 0
                all_data_human_gene_premrna_no_nan.loc[index, on_target_fold] = mean_fold
                all_data_human_gene_premrna_no_nan.loc[index, on_target_fold_normalized] = mean_fold / l
    return all_data_human_gene_premrna_no_nan


def populate_features(df, features, **kwargs):
    if 'self_energy' in features:
        df.loc[:, 'self_energy'] = [float(primer3.calc_homodimer(antisense).dg) for
                                    antisense in
                                    df[SEQUENCE]]
        df.loc[:, 'self_energy'] = df.loc[:,
                                   'self_energy'].astype(float)

    if 'internal_fold' in features:
        df.loc[:, 'internal_fold'] = [RNA.fold(antisense)[1] for antisense in df[SEQUENCE]]

    if 'gc_content' in features:
        df.loc[:, 'gc_content'] = [get_gc_content(seq) for seq in df[SEQUENCE]]

    if 'gc_content_5_prime_5' in features:
        df.loc[:, 'gc_content_5_prime_5'] = [get_gc_content(sequence[-5:]) for sequence in df[SEQUENCE]]

    if 'gc_content_3_prime_5' in features:
        df.loc[:, 'gc_content_3_prime_5'] = [get_gc_content(sequence[:5]) for sequence in df[SEQUENCE]]

    if 'first_nucleotide' in features:
        df.loc[:, 'first_nucleotide'] = [sequence[0] for sequence in df[SEQUENCE]]
    if 'second_nucleotide' in features:
        df.loc[:, 'second_nucleotide'] = [sequence[1] for sequence in df[SEQUENCE]]

    if 'mrna_length' in features:
        df.loc[:, 'mrna_length'] = [len(kwargs['gene_to_data'][gene].full_mrna) for gene in df[CANONICAL_GENE]]
    if 'normalized_start' in features:
        df.loc[:, 'mrna_length'] = [len(kwargs['gene_to_data'][gene].full_mrna) for gene in df[CANONICAL_GENE]]
        df.loc[:, 'normalized_start'] = df[SENSE_START] / df['mrna_length']

    if 'palindromic_fraction' in features:
        df.loc[:, 'palindromic_fraction'] = [palindromic_fraction(seq, 5) for seq in df[SEQUENCE]]

    if 'homooligo_count' in features:
        df.loc[:, 'homooligo_count'] = [homooligo_count(seq) for seq in df[SEQUENCE]]

    if 'entropy' in features:
        df.loc[:, 'entropy'] = [seq_entropy(seq) for seq in df[SEQUENCE]]

    if 'hairpin_score' in features:
        df.loc[:, 'hairpin_score'] = [hairpin_score(seq) for seq in df[SEQUENCE]]

    if 'gc_skew' in features:
        df.loc[:, 'gc_skew'] = [gc_skew(seq) for seq in df[SEQUENCE]]

    if 'at_skew' in features:
        df.loc[:, 'at_skew'] = [at_skew(seq) for seq in df[SEQUENCE]]

    if 'nucleotide_diversity' in features:
        df.loc[:, 'nucleotide_diversity'] = [nucleotide_diversity(seq) for seq in df[SEQUENCE]]

    if 'stop_codon_count' in features:
        df.loc[:, 'stop_codon_count'] = [stop_codon_count(seq) for seq in df[SEQUENCE]]

    if 'at_rich_region_score' in features:
        df.loc[:, 'at_rich_region_score'] = [at_rich_region_score(seq) for seq in df[SEQUENCE]]

    if 'poly_pyrimidine_stretch' in features:
        df.loc[:, 'poly_pyrimidine_stretch'] = [poly_pyrimidine_stretch(seq) for seq in df[SEQUENCE]]
