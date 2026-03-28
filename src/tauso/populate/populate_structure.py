import numpy as np

from ..data.consts import *
from ..features.names import *


def _build_multilength_index(text: bytes, lengths):
    # Iterate backwards so the earliest occurrence overwrites the later ones
    return {L: {text[i : i + L]: i for i in range(len(text) - L, -1, -1)} for L in set(lengths)}


def _find_all_indices(big_string: bytes, small_strings: list[bytes]):
    lengths = [len(s) for s in small_strings]
    indexes = _build_multilength_index(big_string, lengths)
    return np.array([indexes[len(s)].get(s, -1) for s in small_strings], dtype=np.int32)


COMP_U_TABLE = bytes.maketrans(b"ACGTUacgtu", b"UGCAaugcaa")


def get_antisense_u(seq: str) -> str:
    return seq.encode().translate(COMP_U_TABLE)[::-1].decode()


def get_populated_df_with_structure_features(df, genes_u, gene_to_data, use_mask=True):
    if use_mask:
        mask = (df[CELL_LINE_ORGANISM] == "human") & (df[INHIBITION].notna()) & (df[CANONICAL_GENE].isin(genes_u))

        all_data = df[mask].copy()
    else:
        all_data = df.copy()

    n_rows = len(all_data)

    trans_table = str.maketrans("tT", "uU")

    unique_seqs = all_data[SEQUENCE].unique()
    sense_cache = {seq: get_antisense_u(seq) for seq in unique_seqs}

    # OPTIMIZATION 2: Map to a temp column to avoid slow .loc index lookups in the loop
    all_data["__temp_sense"] = all_data[SEQUENCE].map(sense_cache)
    seq_lengths = all_data[SEQUENCE].str.len().values

    pre_mrna_cache = {g: gene_to_data[g].full_mrna.upper().translate(trans_table) for g in genes_u if g in gene_to_data}

    out_start = np.full(n_rows, -1, dtype=np.int32)
    out_start_end = np.zeros(n_rows, dtype=np.int32)
    out_exon = np.zeros(n_rows, dtype=np.int8)
    out_intron = np.zeros(n_rows, dtype=np.int8)
    out_utr = np.zeros(n_rows, dtype=np.int8)
    out_type = np.full(n_rows, "NA", dtype=object)

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

        # Build index, use it, immediately free it
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

        # Original coordinate calculation
        if locus_info.strand == "-":
            gen_coords = locus_info.gene_end - 1 - v_idxs
        else:
            gen_coords = locus_info.gene_start + v_idxs

        coords_2d = gen_coords[:, None]

        exons = np.array(locus_info.exon_indices)
        introns = np.array(locus_info.intron_indices)
        utrs = np.array(locus_info.utr_indices)

        # Original 2D broadcasting
        is_exon = (
            ((coords_2d >= exons[:, 0]) & (coords_2d <= exons[:, 1])).any(axis=1)
            if len(exons) > 0
            else np.zeros(len(gen_coords), dtype=bool)
        )
        is_intron = (
            ((coords_2d >= introns[:, 0]) & (coords_2d <= introns[:, 1])).any(axis=1)
            if len(introns) > 0
            else np.zeros(len(gen_coords), dtype=bool)
        )
        is_utr = (
            ((coords_2d >= utrs[:, 0]) & (coords_2d <= utrs[:, 1])).any(axis=1)
            if len(utrs) > 0
            else np.zeros(len(gen_coords), dtype=bool)
        )

        out_type[v_row_idxs] = "intron"

        utr_mask = is_utr
        exon_mask = is_exon & ~is_utr
        intron_mask = is_intron & ~is_exon

        out_type[v_row_idxs[exon_mask]] = "exon"
        out_type[v_row_idxs[intron_mask]] = "intron"
        out_type[v_row_idxs[utr_mask]] = "utr"

        out_exon[v_row_idxs[exon_mask]] = 1
        out_intron[v_row_idxs[intron_mask]] = 1
        out_utr[v_row_idxs[utr_mask]] = 1

    all_data[SENSE_START] = out_start
    all_data[SENSE_START_FROM_END] = out_start_end
    all_data[SENSE_LENGTH] = seq_lengths
    all_data[SENSE_EXON] = out_exon
    all_data[SENSE_INTRON] = out_intron
    all_data[SENSE_UTR] = out_utr
    all_data[SENSE_TYPE] = out_type

    # Clean up both temp columns
    all_data.drop(columns=["__temp_idx", "__temp_sense"], inplace=True)
    return all_data
