"""Adapt the processed ASOptimizer dataset into OligoAI's input format.

Reads ``PROCESSED_CSV``, maps ASOptimizer columns onto the OligoAI schema, locates
each ASO's sense start in its mRNA, extracts the flanking RNA context, and writes
``ADAPTED_CSV`` for scoring with OligoAI.

Usage:
    python -m notebooks.data.ASOptimizer.to_oligoai_format
"""
import logging

import pandas as pd
from notebooks.data.ASOptimizer.consts import (
    ADAPTED_CSV,
    CELL_LINE_ASOPT,
    CHEMICAL_PATTERN_ASOPT,
    CONTEXT_FLANK,
    INHIBITION_ASOPT,
    LINKAGE_ASOPT,
    LINKAGE_LOCATION_ASOPT,
    PROCESSED_CSV,
    SEQUENCE_ASOPT,
    VOLUME_ASOPT,
)
from notebooks.data.ASOptimizer.utility import transfection_to_oligo
from notebooks.notebook_utils import get_unique_genes

from tauso.common.modifications import transform_linkage_to_oligo, transform_pattern_to_oligo
from tauso.data.consts import (
    BACKBONE_MODS,
    CANONICAL_GENE_NAME,
    INHIBITION_PERCENT,
    STRUCTURE_SENSE_START,
    SUGAR_MODS,
)
from tauso.genome.read_human_genome import get_locus_to_data_dict
from tauso.util import get_antisense_rna

logger = logging.getLogger(__name__)

PHOSPHOROTHIOATE_LINKAGES = ["phosphorothioate", "phosphorothioate/phosphodiester"]


def _mrna(gene, gene_to_data):
    return gene_to_data[gene].full_mrna.upper().replace("T", "U")


def find_sense_start(row, gene_to_data):
    gene = row[CANONICAL_GENE_NAME]
    if pd.isna(gene) or gene not in gene_to_data:
        return -1
    return _mrna(gene, gene_to_data).find(get_antisense_rna(row[SEQUENCE_ASOPT]))


def extract_rna_context(row, gene_to_data):
    mrna = _mrna(row[CANONICAL_GENE_NAME], gene_to_data)
    start = row[STRUCTURE_SENSE_START]
    left = max(0, start - CONTEXT_FLANK)
    right = min(len(mrna), start + len(str(row[SEQUENCE_ASOPT])) + CONTEXT_FLANK)
    return mrna[left:right]


def to_oligoai_format(df, gene_to_data):
    df = df.copy()
    df["split"] = "test"
    df[INHIBITION_PERCENT] = df[INHIBITION_ASOPT]
    df["aso_sequence_5_to_3"] = df[SEQUENCE_ASOPT]
    df["dosage"] = df[VOLUME_ASOPT]
    df[SUGAR_MODS] = df[CHEMICAL_PATTERN_ASOPT].apply(transform_pattern_to_oligo)
    df[BACKBONE_MODS] = df.apply(
        lambda r: transform_linkage_to_oligo(r[LINKAGE_LOCATION_ASOPT], len(r[SEQUENCE_ASOPT])), axis=1
    )

    df[STRUCTURE_SENSE_START] = df.apply(lambda r: find_sense_start(r, gene_to_data), axis=1)
    df = df[df[STRUCTURE_SENSE_START] != -1].copy()

    df["rna_context"] = df.apply(lambda r: extract_rna_context(r, gene_to_data), axis=1)
    df["custom_id"] = df[CANONICAL_GENE_NAME].astype(str) + "_" + df[CELL_LINE_ASOPT].astype(str)
    return transfection_to_oligo(df)


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s | %(message)s")
    df = pd.read_csv(PROCESSED_CSV)
    df = df[df[LINKAGE_ASOPT].isin(PHOSPHOROTHIOATE_LINKAGES)].copy()
    logger.info("loaded %d rows", len(df))

    gene_to_data = get_locus_to_data_dict(include_introns=True, gene_subset=get_unique_genes(df))
    adapted = to_oligoai_format(df, gene_to_data)

    logger.info("writing %d rows -> %s", len(adapted), ADAPTED_CSV)
    adapted.to_csv(ADAPTED_CSV)


if __name__ == "__main__":
    main()
