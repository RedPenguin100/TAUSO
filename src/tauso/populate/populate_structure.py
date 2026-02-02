import numpy as np

from ..data.consts import *
from ..features.names import *
from ..util import get_antisense


def get_populated_df_with_structure_features(df, genes_u, gene_to_data):
    """
    Populate "the data" df with features like exon/intron, start of the sense strand, if found.
    """
    df_copy = df.copy()
    all_data_human = df_copy[df_copy[CELL_LINE_ORGANISM] == 'human']
    all_data_human_no_nan = all_data_human.dropna(subset=[INHIBITION]).copy()
    all_data_human_gene = all_data_human_no_nan[all_data_human_no_nan[CANONICAL_GENE].isin(genes_u)].copy()

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

        idx = pre_mrna.upper().replace("T", "U").find(sense.upper().replace("T", "U"))

        all_data_human_gene.loc[index, SENSE_START] = idx
        all_data_human_gene.loc[index, SENSE_START_FROM_END] = np.abs(
            locus_info.exon_indices[-1][1] - locus_info.gene_start - idx)
        all_data_human_gene.loc[index, SENSE_LENGTH] = len(antisense)
        if idx != -1:
            genome_corrected_index = idx + locus_info.gene_start
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
