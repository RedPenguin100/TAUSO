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
    out_start = np.full(n_rows, -1, dtype=np.int64)
    out_start_end = np.zeros(n_rows, dtype=np.int64)
    out_type = np.full(n_rows, "unannotated", dtype=object)  # Initialize all as unannotated

    out_exon = np.zeros(n_rows, dtype=np.int64)
    out_exon_non_exclusive = np.zeros(n_rows, dtype=np.int64)
    out_intron = np.zeros(n_rows, dtype=np.int64)
    out_3utr = np.zeros(n_rows, dtype=np.int64)
    out_5utr = np.zeros(n_rows, dtype=np.int64)
    out_cds = np.zeros(n_rows, dtype=np.int64)
    out_cds_non_exclusive = np.zeros(n_rows, dtype=np.int64)

    # Normalized sense positions and codon-distance features.
    # NaN for unmapped rows (no gene match), for rows on genes without codon
    # annotations (lncRNAs), and for mRNA-distance features on rows whose
    # genomic coord lies in an intron (no mRNA position exists there).
    out_start_norm = np.full(n_rows, np.nan, dtype=np.float64)
    out_start_end_norm = np.full(n_rows, np.nan, dtype=np.float64)
    out_dist_canonical_stop = np.full(n_rows, np.nan, dtype=np.float64)
    out_dist_closest_stop = np.full(n_rows, np.nan, dtype=np.float64)
    out_dist_canonical_start = np.full(n_rows, np.nan, dtype=np.float64)
    out_dist_closest_start = np.full(n_rows, np.nan, dtype=np.float64)
    out_mrna_dist_canonical_stop = np.full(n_rows, np.nan, dtype=np.float64)
    out_mrna_dist_closest_stop = np.full(n_rows, np.nan, dtype=np.float64)

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

        # Normalized positions: sense_start / gene_length, sense_start_from_end / gene_length.
        # Range [0, 1]; tells the model "how far into the pre-mRNA is the target?"
        # length-invariant across genes (PSD3 ~640k vs APOL1 ~15k).
        gene_length = float(locus_info.gene_end - locus_info.gene_start)
        if gene_length > 0:
            out_start_norm[v_row_idxs] = v_idxs.astype(np.float64) / gene_length
            out_start_end_norm[v_row_idxs] = (gene_length - v_idxs).astype(np.float64) / gene_length

        # Genomic distance from the ASO center to the nearest canonical /
        # any-transcript stop or start codon. NaN for genes with no codon
        # annotation (lncRNAs etc). Distance is in genomic bp (intron-inclusive);
        # mRNA-distance variants live below.
        if locus_info.stop_codons:
            stop_mids = np.array([(s + e) / 2.0 for s, e in locus_info.stop_codons], dtype=np.float64)
            out_dist_canonical_stop[v_row_idxs] = np.min(
                np.abs(gen_coords[:, None] - stop_mids[None, :]), axis=1
            )
        if locus_info.all_stop_codons:
            all_stop_mids = np.array([(s + e) / 2.0 for s, e in locus_info.all_stop_codons], dtype=np.float64)
            out_dist_closest_stop[v_row_idxs] = np.min(
                np.abs(gen_coords[:, None] - all_stop_mids[None, :]), axis=1
            )
        if locus_info.start_codons:
            start_mids = np.array([(s + e) / 2.0 for s, e in locus_info.start_codons], dtype=np.float64)
            out_dist_canonical_start[v_row_idxs] = np.min(
                np.abs(gen_coords[:, None] - start_mids[None, :]), axis=1
            )
        if locus_info.all_start_codons:
            all_start_mids = np.array([(s + e) / 2.0 for s, e in locus_info.all_start_codons], dtype=np.float64)
            out_dist_closest_start[v_row_idxs] = np.min(
                np.abs(gen_coords[:, None] - all_start_mids[None, :]), axis=1
            )

        # mRNA-distance variants for stop codons. Walk the canonical exon set in
        # 5'->3' order and map a genomic coord to its position in the spliced
        # mRNA. Intronic positions return NaN (no mRNA position exists).
        if locus_info.stop_codons or locus_info.all_stop_codons:
            sorted_exons = sorted(
                locus_info._exon_indices,
                key=(lambda x: -x[0]) if locus_info.strand == StrandType.NEG else (lambda x: x[0]),
            )
            cum_pre = np.zeros(len(sorted_exons), dtype=np.int64)
            running = 0
            for i, (s, e) in enumerate(sorted_exons):
                cum_pre[i] = running
                running += e - s

            def _gen_to_mrna(coords, sorted_exons=sorted_exons, cum_pre=cum_pre, locus_info=locus_info):
                """Vectorized: genomic coords → mRNA positions (NaN if intronic)."""
                out_mrna = np.full(len(coords), np.nan, dtype=np.float64)
                for (s, e), pre in zip(sorted_exons, cum_pre):
                    mask = (coords >= s) & (coords < e)
                    if not np.any(mask):
                        continue
                    if locus_info.strand == StrandType.NEG:
                        out_mrna[mask] = pre + (e - 1 - coords[mask])
                    else:
                        out_mrna[mask] = pre + (coords[mask] - s)
                return out_mrna

            aso_mrna = _gen_to_mrna(gen_coords.astype(np.float64))

            if locus_info.stop_codons:
                stop_mids_arr = np.array([(s + e) / 2.0 for s, e in locus_info.stop_codons], dtype=np.float64)
                stop_mrna = _gen_to_mrna(stop_mids_arr)
                # min over stops, NaN-aware
                valid_stops = stop_mrna[~np.isnan(stop_mrna)]
                if len(valid_stops) > 0:
                    dists = np.min(np.abs(aso_mrna[:, None] - valid_stops[None, :]), axis=1)
                    out_mrna_dist_canonical_stop[v_row_idxs] = dists
            if locus_info.all_stop_codons:
                all_stop_mids_arr = np.array([(s + e) / 2.0 for s, e in locus_info.all_stop_codons], dtype=np.float64)
                all_stop_mrna = _gen_to_mrna(all_stop_mids_arr)
                valid_all_stops = all_stop_mrna[~np.isnan(all_stop_mrna)]
                if len(valid_all_stops) > 0:
                    dists_all = np.min(np.abs(aso_mrna[:, None] - valid_all_stops[None, :]), axis=1)
                    out_mrna_dist_closest_stop[v_row_idxs] = dists_all

        # Much cleaner interval checking
        is_exon = _in_intervals(gen_coords, locus_info._exon_indices)
        is_intron = _in_intervals(gen_coords, locus_info._intron_indices)
        is_3utr = _in_intervals(gen_coords, locus_info._3utr_indices)
        is_5utr = _in_intervals(gen_coords, locus_info._5utr_indices)
        is_generic_utr = _in_intervals(gen_coords, locus_info.utr_indices)

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
        out_exon_non_exclusive[v_row_idxs[is_exon]] = 1
        out_intron[v_row_idxs[mask_intron]] = 1
        out_3utr[v_row_idxs[mask_3utr]] = 1
        out_5utr[v_row_idxs[mask_5utr]] = 1

        # sense_cds and sense_cds_non_exclusive = exon hit AND the gene is
        # protein-coding, in exclusive (mask_exon, UTR excluded) vs non-exclusive
        # (is_exon, any exonic hit including UTRs that overlap exons) forms.
        # The biotype guard keeps both at zero on lncRNAs (MALAT1, SNHG14, …)
        # where no UTRs are annotated to mask and every exonic hit would
        # otherwise spuriously look like CDS.
        if locus_info.gene_type == GeneType.PROTEIN_CODING:
            out_cds[v_row_idxs[mask_exon]] = 1
            out_cds_non_exclusive[v_row_idxs[is_exon]] = 1

    # Apply exclusive features
    all_data[SENSE_START] = out_start
    all_data[SENSE_START_FROM_END] = out_start_end
    all_data[SENSE_LENGTH] = seq_lengths
    all_data[SENSE_EXON] = out_exon
    all_data[SENSE_EXON_NON_EXCLUSIVE] = out_exon_non_exclusive
    all_data[SENSE_INTRON] = out_intron
    all_data[SENSE_UTR] = out_3utr | out_5utr
    all_data[SENSE_3UTR] = out_3utr
    all_data[SENSE_5UTR] = out_5utr
    all_data[SENSE_CDS] = out_cds
    all_data[SENSE_CDS_NON_EXCLUSIVE] = out_cds_non_exclusive
    all_data[SENSE_TYPE] = out_type

    # Normalized positions + codon distances (genomic + mRNA)
    all_data[SENSE_START_NORM] = out_start_norm
    all_data[SENSE_START_FROM_END_NORM] = out_start_end_norm
    all_data[SENSE_DIST_TO_CANONICAL_STOP] = out_dist_canonical_stop
    all_data[SENSE_DIST_TO_CLOSEST_STOP] = out_dist_closest_stop
    all_data[SENSE_DIST_TO_CANONICAL_START] = out_dist_canonical_start
    all_data[SENSE_DIST_TO_CLOSEST_START] = out_dist_closest_start
    all_data[SENSE_MRNA_DIST_TO_CANONICAL_STOP] = out_mrna_dist_canonical_stop
    all_data[SENSE_MRNA_DIST_TO_CLOSEST_STOP] = out_mrna_dist_closest_stop

    all_data.drop(columns=["__temp_idx", "__temp_sense"], inplace=True)
    return all_data
