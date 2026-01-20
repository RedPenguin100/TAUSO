from notebooks.preprocessing import get_unique_genes
from tauso.new_model.consts_dataframe import CELL_LINE_ORGANISM, INHIBITION, CANONICAL_GENE

import pytest

import pandas as pd
import numpy as np

from tauso.consts import DATA_PATH

@pytest.mark.skip
def test_genes_unique():
    all_data = pd.read_csv(DATA_PATH / 'data_from_article_fixed.csv')
    genes_u = get_unique_genes(all_data)

    assert 'KRAS' in genes_u
    assert 'HBV' not in genes_u
    assert len(set(genes_u)) == len(genes_u)

# TODO: improve test
# def test_intron_annotation():
#     all_data = pd.read_csv(DATA_PATH / 'data_from_article_fixed.csv')
#     all_data_no_nan = all_data.dropna(subset=[INHIBITION]).copy()
#
#     genes_u = get_unique_human_genes(all_data)
#
#     all_data_no_nan_human = all_data_no_nan[all_data_no_nan[CELL_LINE_ORGANISM] == 'human']
#     all_data_human_gene = all_data_no_nan_human[all_data_no_nan_human[CANONICAL_GENE].isin(genes_u)].copy()
#
#     gene_to_data = get_gene_to_data(genes_u)
#
#     found = False
#
#     all_data_human_gene[SENSE_START] = np.zeros_like(all_data_human_gene[CANONICAL_GENE], dtype=int)
#     all_data_human_gene[SENSE_LENGTH] = np.zeros_like(all_data_human_gene[CANONICAL_GENE], dtype=int)
#     all_data_human_gene[SENSE_TYPE] = "NA"
#     for index, row in all_data_human_gene.iterrows():
#         gene_name = row[CANONICAL_GENE]
#         locus_info = gene_to_data[gene_name]
#         pre_mrna = locus_info.full_mrna
#         antisense = row[SEQUENCE]
#         sense = get_antisense(antisense)
#         idx = pre_mrna.find(sense)
#         all_data_human_gene.loc[index, SENSE_START] = idx
#         all_data_human_gene.loc[index, SENSE_LENGTH] = len(antisense)
#         if idx != -1:
#             genome_corrected_index = idx + locus_info.exon_indices[0][0]
#             found = False
#             for exon_indices in locus_info.exon_indices:
#                 # print(exon[0], exon[1])
#                 if exon_indices[0] <= genome_corrected_index <= exon_indices[1]:
#                     all_data_human_gene.loc[index, SENSE_TYPE] = 'exon'
#                     found = True
#                     break
#         if not found:
#             all_data_human_gene.loc[index, SENSE_TYPE] = 'intron'
