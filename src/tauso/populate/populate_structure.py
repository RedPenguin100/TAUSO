import logging

import numpy as np

from ..data.consts import *
from ..genome.LocusInfo import GeneType, StrandType
from ..util import get_antisense_rna

logger = logging.getLogger(__name__)


def _build_multilength_index(text: bytes, lengths: set):
    """Builds a lookup index for finding sequence matches efficiently."""
    return {L: {text[i : i + L]: i for i in range(len(text) - L, -1, -1)} for L in lengths}


def _in_intervals(coords: np.ndarray, intervals: list) -> np.ndarray:
    """Returns a boolean array indicating if coords fall within any [start, end) interval."""
    if not intervals or len(intervals) == 0:
        return np.zeros(len(coords), dtype=bool)

    arr = np.array(intervals)
    return ((coords[:, None] >= arr[:, 0]) & (coords[:, None] < arr[:, 1])).any(axis=1)


def get_populated_df_with_structure_features(df, genes_u, gene_to_data, use_mask=True):
    if use_mask:
        mask = (df[CELL_LINE_ORGANISM] == "human") & (df[CANONICAL_GENE].isin(genes_u))
        all_data = df[mask].copy()
    else:
        all_data = df.copy()

    n_rows = len(all_data)
    if n_rows == 0:
        return all_data

    trans_table = str.maketrans("tT", "uU")

    # Encode sequences to bytes outside the loop for speed
    unique_seqs = all_data[SEQUENCE].unique()
    sense_cache = {seq: get_antisense_rna(seq).encode() for seq in unique_seqs}

    all_data["__temp_sense"] = all_data[SEQUENCE].map(sense_cache)
    seq_lengths = all_data[SEQUENCE].str.len().values

    pre_mrna_cache = {
        g: gene_to_data[g].full_mrna.upper().translate(trans_table).encode() for g in genes_u if g in gene_to_data
    }

    # Initialize tracking arrays
    out_start = np.full(n_rows, -1, dtype=np.int32)
    out_start_end = np.zeros(n_rows, dtype=np.int32)
    out_type = np.full(n_rows, "unannotated", dtype=object)  # Initialize all as unannotated

    out_exon = np.zeros(n_rows, dtype=np.int8)
    out_intron = np.zeros(n_rows, dtype=np.int8)
    out_3utr = np.zeros(n_rows, dtype=np.int8)
    out_5utr = np.zeros(n_rows, dtype=np.int8)
    out_cds = np.zeros(n_rows, dtype=np.int8)

    out_n_exon = np.zeros(n_rows, dtype=np.int8)
    out_n_intron = np.zeros(n_rows, dtype=np.int8)
    out_n_3utr = np.zeros(n_rows, dtype=np.int8)
    out_n_5utr = np.zeros(n_rows, dtype=np.int8)
    out_n_utr = np.zeros(n_rows, dtype=np.int8)

    all_data["__temp_idx"] = np.arange(n_rows)

    for gene_name, group in all_data.groupby(CANONICAL_GENE, observed=True):
        if gene_name not in pre_mrna_cache:
            continue

        locus_info = gene_to_data[gene_name]
        pre_mrna = pre_mrna_cache[gene_name]

        row_idxs = group["__temp_idx"].values
        group_senses_b = group["__temp_sense"].values

        group_lengths = set(len(s) for s in group_senses_b)
        index = _build_multilength_index(pre_mrna, group_lengths)
        idxs = np.array([index[len(s)].get(s, -1) for s in group_senses_b], dtype=np.int32)

        valid_mask = idxs != -1
        if not valid_mask.any():
            continue

        v_row_idxs = row_idxs[valid_mask]
        v_idxs = idxs[valid_mask]

        out_start[v_row_idxs] = v_idxs
        out_start_end[v_row_idxs] = len(pre_mrna) - v_idxs

        # Get lengths of the valid matches
        v_seq_lengths = seq_lengths[v_row_idxs]

        # Calculate the 5'-tilted center offset
        # Len 20 -> offset 9, Len 21 -> offset 10
        center_offsets = (v_seq_lengths - 1) // 2

        if locus_info.strand == StrandType.NEG:
            # Map the center to the negative strand's absolute genomic coordinate
            gen_coords = locus_info.gene_end - 1 - v_idxs - center_offsets
        else:
            # Map the center to the positive strand's absolute genomic coordinate
            gen_coords = locus_info.gene_start + v_idxs + center_offsets

        # Much cleaner interval checking
        is_exon = _in_intervals(gen_coords, locus_info._exon_indices)
        is_intron = _in_intervals(gen_coords, locus_info._intron_indices)
        is_3utr = _in_intervals(gen_coords, locus_info._3utr_indices)
        is_5utr = _in_intervals(gen_coords, locus_info._5utr_indices)
        is_generic_utr = _in_intervals(gen_coords, locus_info.utr_indices)

        # --- NON-EXCLUSIVE ASSIGNMENTS ---
        out_n_exon[v_row_idxs[is_exon]] = 1
        out_n_intron[v_row_idxs[is_intron]] = 1
        out_n_3utr[v_row_idxs[is_3utr]] = 1
        out_n_5utr[v_row_idxs[is_5utr]] = 1
        out_n_utr[v_row_idxs[is_generic_utr | is_3utr | is_5utr]] = 1

        # --- EXCLUSIVE ASSIGNMENTS & TYPE ---
        # Note: generic_utr_mask was completely unused for exclusive flagging, so it's removed.
        mask_3utr = is_3utr
        mask_5utr = is_5utr & ~mask_3utr
        mask_exon = is_exon & ~(mask_3utr | mask_5utr | is_generic_utr)
        mask_intron = is_intron & ~is_exon

        out_type[v_row_idxs[mask_exon]] = "exon"
        out_type[v_row_idxs[mask_intron]] = "intron"
        out_type[v_row_idxs[mask_3utr]] = "3UTR"
        out_type[v_row_idxs[mask_5utr]] = "5UTR"

        out_exon[v_row_idxs[mask_exon]] = 1
        out_intron[v_row_idxs[mask_intron]] = 1
        out_3utr[v_row_idxs[mask_3utr]] = 1
        out_5utr[v_row_idxs[mask_5utr]] = 1

        # sense_cds = sense_exon AND the gene is protein-coding. mask_exon already
        # excludes UTRs, so for protein-coding genes it's the CDS region — but
        # for lncRNAs (MALAT1, SNHG14, ...) every exonic hit also satisfies
        # mask_exon (no UTRs to mask), which would spuriously label them CDS.
        # The biotype guard ensures sense_cds is 1 only for true CDS hits.
        if locus_info.gene_type == GeneType.PROTEIN_CODING:
            out_cds[v_row_idxs[mask_exon]] = 1

    # Apply exclusive features
    all_data[SENSE_START] = out_start
    all_data[SENSE_START_FROM_END] = out_start_end
    all_data[SENSE_LENGTH] = seq_lengths
    all_data[SENSE_EXON] = out_exon
    all_data[SENSE_INTRON] = out_intron
    all_data[SENSE_UTR] = out_3utr | out_5utr
    all_data[SENSE_3UTR] = out_3utr
    all_data[SENSE_5UTR] = out_5utr
    all_data[SENSE_CDS] = out_cds
    all_data[SENSE_TYPE] = out_type

    # Apply non-exclusive features
    all_data[f"n_{SENSE_EXON}"] = out_n_exon
    all_data[f"n_{SENSE_INTRON}"] = out_n_intron
    all_data[f"n_{SENSE_3UTR}"] = out_n_3utr
    all_data[f"n_{SENSE_5UTR}"] = out_n_5utr
    all_data[f"n_{SENSE_UTR}"] = out_n_utr

    all_data.drop(columns=["__temp_idx", "__temp_sense"], inplace=True)
    return all_data
