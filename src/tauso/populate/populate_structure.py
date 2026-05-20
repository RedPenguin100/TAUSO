import logging

import numpy as np

from ..data.consts import *
from ..features.names import *
from ..genome.LocusInfo import StrandType
from ..util import get_antisense_u

logger = logging.getLogger(__name__)


def _build_multilength_index(text: bytes, lengths):
    # Iterate backwards so the earliest occurrence overwrites the later ones
    return {L: {text[i : i + L]: i for i in range(len(text) - L, -1, -1)} for L in set(lengths)}


def _find_all_indices(big_string: bytes, small_strings: list[bytes]):
    lengths = [len(s) for s in small_strings]
    indexes = _build_multilength_index(big_string, lengths)
    return np.array([indexes[len(s)].get(s, -1) for s in small_strings], dtype=np.int32)


def get_populated_df_with_structure_features(df, genes_u, gene_to_data, use_mask=True):
    if use_mask:
        mask = (df[CELL_LINE_ORGANISM] == "human") & (df[CANONICAL_GENE].isin(genes_u))
        all_data = df[mask].copy()
    else:
        all_data = df.copy()

    n_rows = len(all_data)
    trans_table = str.maketrans("tT", "uU")

    unique_seqs = all_data[SEQUENCE].unique()
    sense_cache = {seq: get_antisense_u(seq) for seq in unique_seqs}

    all_data["__temp_sense"] = all_data[SEQUENCE].map(sense_cache)
    seq_lengths = all_data[SEQUENCE].str.len().values

    pre_mrna_cache = {g: gene_to_data[g].full_mrna.upper().translate(trans_table) for g in genes_u if g in gene_to_data}

    # Initialize exclusive arrays
    out_start = np.full(n_rows, -1, dtype=np.int32)
    out_start_end = np.zeros(n_rows, dtype=np.int32)
    out_exon = np.zeros(n_rows, dtype=np.int8)
    out_intron = np.zeros(n_rows, dtype=np.int8)
    out_utr3 = np.zeros(n_rows, dtype=np.int8)
    out_utr5 = np.zeros(n_rows, dtype=np.int8)
    out_type = np.full(n_rows, "NA", dtype=object)

    # Initialize NON-exclusive arrays
    out_n_exon = np.zeros(n_rows, dtype=np.int8)
    out_n_intron = np.zeros(n_rows, dtype=np.int8)
    out_n_utr3 = np.zeros(n_rows, dtype=np.int8)
    out_n_utr5 = np.zeros(n_rows, dtype=np.int8)
    out_n_utr = np.zeros(n_rows, dtype=np.int8)

    all_data["__temp_idx"] = np.arange(n_rows)

    for gene_name, group in all_data.groupby(CANONICAL_GENE):
        if gene_name not in pre_mrna_cache:
            continue

        locus_info = gene_to_data[gene_name]
        pre_mrna_u = pre_mrna_cache[gene_name]
        pre_mrna_len = len(pre_mrna_u)

        row_idxs = group["__temp_idx"].values
        group_senses = group["__temp_sense"].values
        group_senses_b = [s.encode() for s in group_senses]

        group_lengths = set(len(s) for s in group_senses_b)
        index = _build_multilength_index(pre_mrna_u.encode(), group_lengths)
        idxs = np.array([index[len(s)].get(s, -1) for s in group_senses_b], dtype=np.int32)

        valid_mask = idxs != -1
        if not valid_mask.any():
            continue

        v_row_idxs = row_idxs[valid_mask]
        v_idxs = idxs[valid_mask]

        out_start[v_row_idxs] = v_idxs
        out_start_end[v_row_idxs] = pre_mrna_len - v_idxs

        if locus_info.strand == StrandType.NEG:
            gen_coords = locus_info.gene_end - 1 - v_idxs
        else:
            gen_coords = locus_info.gene_start + v_idxs

        coords_2d = gen_coords[:, None]

        exons = np.array(locus_info._exon_indices)
        introns = np.array(locus_info._intron_indices)
        utrs3 = np.array(locus_info._3utr_indices)
        utrs5 = np.array(locus_info._5utr_indices)
        utrs = np.array(locus_info.utr_indices)  # For capturing all UTRs

        is_exon = (
            ((coords_2d >= exons[:, 0]) & (coords_2d < exons[:, 1])).any(axis=1)
            if len(exons) > 0
            else np.zeros(len(gen_coords), dtype=bool)
        )
        is_intron = (
            ((coords_2d >= introns[:, 0]) & (coords_2d < introns[:, 1])).any(axis=1)
            if len(introns) > 0
            else np.zeros(len(gen_coords), dtype=bool)
        )
        is_utr3 = (
            ((coords_2d >= utrs3[:, 0]) & (coords_2d < utrs3[:, 1])).any(axis=1)
            if len(utrs3) > 0
            else np.zeros(len(gen_coords), dtype=bool)
        )
        is_utr5 = (
            ((coords_2d >= utrs5[:, 0]) & (coords_2d < utrs5[:, 1])).any(axis=1)
            if len(utrs5) > 0
            else np.zeros(len(gen_coords), dtype=bool)
        )
        is_generic_utr = (
            ((coords_2d >= utrs[:, 0]) & (coords_2d < utrs[:, 1])).any(axis=1)
            if len(utrs) > 0
            else np.zeros(len(gen_coords), dtype=bool)
        )

        # --- NON-EXCLUSIVE ASSIGNMENTS ---
        out_n_exon[v_row_idxs[is_exon]] = 1
        out_n_intron[v_row_idxs[is_intron]] = 1
        out_n_utr3[v_row_idxs[is_utr3]] = 1
        out_n_utr5[v_row_idxs[is_utr5]] = 1
        out_n_utr[v_row_idxs[is_generic_utr | is_utr3 | is_utr5]] = 1

        # --- EXCLUSIVE ASSIGNMENTS & TYPE ---
        out_type[v_row_idxs] = "unannotated"

        utr3_mask = is_utr3
        utr5_mask = is_utr5 & ~is_utr3
        generic_utr_mask = is_generic_utr & ~(is_utr3 | is_utr5)
        exon_mask = is_exon & ~(is_utr3 | is_utr5 | is_generic_utr)
        intron_mask = is_intron & ~is_exon

        out_type[v_row_idxs[exon_mask]] = "exon"
        out_type[v_row_idxs[intron_mask]] = "intron"
        out_type[v_row_idxs[utr3_mask]] = "3UTR"
        out_type[v_row_idxs[utr5_mask]] = "5UTR"

        out_exon[v_row_idxs[exon_mask]] = 1
        out_intron[v_row_idxs[intron_mask]] = 1
        out_utr3[v_row_idxs[utr3_mask]] = 1
        out_utr5[v_row_idxs[utr5_mask]] = 1

    # Apply exclusive features
    all_data[SENSE_START] = out_start
    all_data[SENSE_START_FROM_END] = out_start_end
    all_data[SENSE_LENGTH] = seq_lengths
    all_data[SENSE_EXON] = out_exon
    all_data[SENSE_INTRON] = out_intron
    all_data[SENSE_UTR] = out_utr3 | out_utr5
    all_data[SENSE_3UTR] = out_utr3
    all_data[SENSE_5UTR] = out_utr5
    all_data[SENSE_TYPE] = out_type

    # Apply non-exclusive features
    all_data[f"n_{SENSE_EXON}"] = out_n_exon
    all_data[f"n_{SENSE_INTRON}"] = out_n_intron
    all_data[f"n_{SENSE_3UTR}"] = out_n_utr3
    all_data[f"n_{SENSE_5UTR}"] = out_n_utr5
    all_data[f"n_{SENSE_UTR}"] = out_n_utr

    all_data.drop(columns=["__temp_idx", "__temp_sense"], inplace=True)
    return all_data
