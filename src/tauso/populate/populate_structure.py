import logging
from types import SimpleNamespace

import numpy as np

from ..data.consts import *
from ..genome.LocusInfo import GeneType, StrandType
from ..util import get_antisense_rna

logger = logging.getLogger(__name__)


def _build_multilength_index(text: bytes, lengths: set):
    """text: bytes, lengths: set of int -> {L: {substring: position}} (last occurrence wins), per L."""
    return {L: {text[i : i + L]: i for i in range(len(text) - L, -1, -1)} for L in lengths}


def _in_intervals(coords: np.ndarray, intervals: list) -> np.ndarray:
    """coords: int array, intervals: list of (start, end) -> bool array: is each coord in any interval."""
    inside = np.zeros(len(coords), dtype=bool)
    for s, e in intervals:
        inside |= (coords >= s) & (coords < e)
    return inside


def _host_interval_len(coords: np.ndarray, intervals: list) -> np.ndarray:
    """coords: int array, intervals: list of (start, end) -> float array: length of the interval each coord
    is in, else 0."""
    length = np.zeros(len(coords), dtype=np.float64)
    for s, e in intervals:
        length[(coords >= s) & (coords < e)] = e - s
    return length


def _codon_midpoints(codons) -> np.ndarray:
    """codons: list of (start, end) -> float array of interval midpoints."""
    return np.array([(s + e) / 2.0 for s, e in codons], dtype=np.float64)


def _genomic_to_mrna(coords, sorted_exons, cum_pre, neg_strand):
    """Map genomic coords to spliced-mRNA positions; NaN for intronic coords."""
    out = np.full(len(coords), np.nan, dtype=np.float64)
    for (s, e), pre in zip(sorted_exons, cum_pre):
        mask = (coords >= s) & (coords < e)
        out[mask] = pre + (e - 1 - coords[mask] if neg_strand else coords[mask] - s)
    return out


def _validate_gene_lengths(genes, gene_to_data):
    """Fail loudly if any processed gene has non-positive genomic length (would break normalization)."""
    bad = [g for g in genes if gene_to_data[g].gene_end <= gene_to_data[g].gene_start]
    if bad:
        raise ValueError(f"genes with non-positive genomic length: {bad}")


def _validate_codons(genes, gene_to_data):
    """Fail loudly if any processed protein-coding gene lacks a canonical start or stop codon (broken CDS
    annotation -- such a gene should not be an ASO target)."""
    bad = [
        g
        for g in genes
        if gene_to_data[g].gene_type == GeneType.PROTEIN_CODING
        and (gene_to_data[g].stop_codon is None or gene_to_data[g].start_codon is None)
    ]
    if bad:
        raise ValueError(f"protein-coding genes missing a canonical start/stop codon (broken annotation): {bad}")


def _match_positions(pre_mrna, senses_b) -> np.ndarray:
    """Start position of each antisense-RNA sequence within the pre-mRNA (bytes), or -1 if absent."""
    index = _build_multilength_index(pre_mrna, {len(s) for s in senses_b})
    return np.array([index[len(s)].get(s, -1) for s in senses_b], dtype=np.int32)


def _genetic_coordinates(hit_pos, match_lengths, locus_info) -> np.ndarray:
    """Genomic coordinate of each match's 5'-tilted center (len 20 -> +9, len 21 -> +10)."""
    center = (match_lengths - 1) // 2
    if locus_info.strand == StrandType.NEG:
        return locus_info.gene_end - 1 - hit_pos - center
    return locus_info.gene_start + hit_pos + center


# --- per-gene feature writers: each computes one feature block for the matched rows (`hit_rows`) and
#     writes it into the shared outputs object `out`. Rows they don't touch keep their default (NaN / 0). ---


def assign_target_position(out, hit_rows, hit_pos, pre_mrna_len, locus_info):
    """Target start within the pre-mRNA (from each end), raw and normalized to [0, 1] (length-invariant)."""
    gene_length = float(locus_info.gene_end - locus_info.gene_start)
    out.start[hit_rows] = hit_pos
    out.start_end[hit_rows] = pre_mrna_len - hit_pos
    out.start_norm[hit_rows] = hit_pos / gene_length
    out.start_end_norm[hit_rows] = (gene_length - hit_pos) / gene_length


def assign_canonical_splice_junction_distances(out, hit_rows, genetic_coordinates, locus_info):
    """Distance (bp) to the nearest CANONICAL-transcript intron boundary, filed by exon/intron side so a 0
    is unambiguous, both raw and log1p (the model-facing form). Single-exon genes stay NaN."""
    splice_junction = np.array(locus_info._intron_indices).flatten()
    if splice_junction.size == 0:
        return
    dist = np.abs(genetic_coordinates[:, None] - splice_junction).min(axis=1)
    in_exon = _in_intervals(genetic_coordinates, locus_info._exon_indices)
    in_intron = _in_intervals(genetic_coordinates, locus_info._intron_indices)
    out.dist_sj_exonic[hit_rows] = np.where(in_exon, dist, np.nan)
    out.dist_sj_intronic[hit_rows] = np.where(in_exon, np.nan, dist)
    out.junction_logdist_exonic[hit_rows] = np.where(in_exon, np.log1p(dist), np.nan)
    out.junction_logdist_intronic[hit_rows] = np.where(in_intron, np.log1p(dist), np.nan)


def assign_closest_splice_junction_distance(out, hit_rows, genetic_coordinates, locus_info):
    """Distance to the nearest splice junction across ALL transcript isoforms (side-agnostic, single value),
    raw and log1p (the model-facing form). NaN if the gene is single-exon in every isoform."""
    sj = np.asarray(locus_info.all_splice_junctions)
    if sj.size == 0:
        return
    dist = np.abs(genetic_coordinates[:, None] - sj).min(axis=1)
    out.dist_closest_sj[hit_rows] = dist
    out.junction_logdist_closest[hit_rows] = np.log1p(dist)


def assign_host_lengths(out, hit_rows, genetic_coordinates, locus_info):
    """log length of the host exon / intron the coord sits in (exon wins overlaps). NaN if in neither."""
    exon_len = _host_interval_len(genetic_coordinates, locus_info._exon_indices)
    intron_len = _host_interval_len(genetic_coordinates, locus_info._intron_indices)
    in_exon = exon_len > 0
    in_intron = (intron_len > 0) & ~in_exon
    out.host_exon_log_len[hit_rows[in_exon]] = np.log(np.maximum(exon_len[in_exon], 1))
    out.host_intron_log_len[hit_rows[in_intron]] = np.log(np.maximum(intron_len[in_intron], 1))


def assign_distance_to_special_codons(out, hit_rows, genetic_coordinates, locus_info):
    """Genomic (intron-inclusive) distance to the nearest special codon -- start / stop, canonical
    transcript (a single codon) and across all isoforms. Non-coding genes (no such codons) stay NaN."""
    stop, start = locus_info.stop_codon, locus_info.start_codon
    for codons, target in (
        ([stop] if stop else [], out.dist_canonical_stop),
        (locus_info.all_stop_codons, out.dist_closest_stop),
        ([start] if start else [], out.dist_canonical_start),
        (locus_info.all_start_codons, out.dist_closest_start),
    ):
        if codons:
            mids = _codon_midpoints(codons)
            target[hit_rows] = np.min(np.abs(genetic_coordinates[:, None] - mids[None, :]), axis=1)


def assign_mrna_stop_distance(out, hit_rows, genetic_coordinates, locus_info):
    """Spliced-mRNA distance to the nearest stop codon (canonical + all-isoform), walking canonical exons
    5'->3' so introns are excluded. Non-coding genes stay NaN."""
    neg = locus_info.strand == StrandType.NEG
    exons = sorted(locus_info._exon_indices, key=(lambda x: -x[0]) if neg else (lambda x: x[0]))
    cum_pre = np.cumsum([0] + [e - s for s, e in exons])[:-1]  # cumulative exon length before each exon
    aso_mrna = None
    stop = locus_info.stop_codon
    for codons, target in (
        ([stop] if stop else [], out.mrna_dist_canonical_stop),
        (locus_info.all_stop_codons, out.mrna_dist_closest_stop),
    ):
        if not codons:
            continue
        if aso_mrna is None:
            aso_mrna = _genomic_to_mrna(genetic_coordinates.astype(np.float64), exons, cum_pre, neg)
        stop_mrna = _genomic_to_mrna(_codon_midpoints(codons), exons, cum_pre, neg)
        stop_mrna = stop_mrna[~np.isnan(stop_mrna)]
        if len(stop_mrna):
            target[hit_rows] = np.min(np.abs(aso_mrna[:, None] - stop_mrna[None, :]), axis=1)


def assign_region(out, hit_rows, genetic_coordinates, locus_info):
    """Exclusive region assignment: type + in_exon / intron / 3UTR / 5UTR / CDS flags. Priority
    3UTR > 5UTR > exon; intron only outside exons; CDS only on protein-coding genes (keeps lncRNAs at 0)."""
    in_exon = _in_intervals(genetic_coordinates, locus_info._exon_indices)
    in_intron = _in_intervals(genetic_coordinates, locus_info._intron_indices)
    in_3utr = _in_intervals(genetic_coordinates, locus_info._3utr_indices)
    in_5utr = _in_intervals(genetic_coordinates, locus_info._5utr_indices)
    in_generic_utr = _in_intervals(genetic_coordinates, locus_info.utr_indices)

    utr3 = in_3utr
    utr5 = in_5utr & ~utr3
    exon = in_exon & ~(utr3 | utr5 | in_generic_utr)
    intron = in_intron & ~in_exon

    out.type[hit_rows[exon]] = "exon"
    out.type[hit_rows[intron]] = "intron"
    out.type[hit_rows[utr3]] = "3UTR"
    out.type[hit_rows[utr5]] = "5UTR"
    out.exon[hit_rows[exon]] = 1
    out.exon_non_exclusive[hit_rows[in_exon]] = 1
    out.intron[hit_rows[intron]] = 1
    out.utr3[hit_rows[utr3]] = 1
    out.utr5[hit_rows[utr5]] = 1
    if locus_info.gene_type == GeneType.PROTEIN_CODING:
        out.cds[hit_rows[exon]] = 1
        out.cds_non_exclusive[hit_rows[in_exon]] = 1


def _init_outputs(n_rows):
    """One array per feature. Float features default to NaN (unmapped rows, genes with no codons, ...);
    int flags default to 0; start defaults to -1."""

    def nan():
        return np.full(n_rows, np.nan, dtype=np.float64)

    return SimpleNamespace(
        start=np.full(n_rows, -1, dtype=np.int64),
        start_end=np.zeros(n_rows, dtype=np.int64),
        type=np.full(n_rows, "unannotated", dtype=object),
        exon=np.zeros(n_rows, dtype=np.int8),
        exon_non_exclusive=np.zeros(n_rows, dtype=np.int8),
        intron=np.zeros(n_rows, dtype=np.int8),
        utr3=np.zeros(n_rows, dtype=np.int8),
        utr5=np.zeros(n_rows, dtype=np.int8),
        cds=np.zeros(n_rows, dtype=np.int8),
        cds_non_exclusive=np.zeros(n_rows, dtype=np.int8),
        start_norm=nan(),
        start_end_norm=nan(),
        dist_canonical_stop=nan(),
        dist_closest_stop=nan(),
        dist_canonical_start=nan(),
        dist_closest_start=nan(),
        mrna_dist_canonical_stop=nan(),
        mrna_dist_closest_stop=nan(),
        dist_sj_exonic=nan(),
        dist_sj_intronic=nan(),
        dist_closest_sj=nan(),
        host_exon_log_len=nan(),
        host_intron_log_len=nan(),
        junction_logdist_exonic=nan(),
        junction_logdist_intronic=nan(),
        junction_logdist_closest=nan(),
    )


def get_populated_df_with_structure_features(df, genes_u, gene_to_data, use_mask=True):
    if use_mask:
        mask = (df[CELL_LINE_ORGANISM] == "human") & (df[CANONICAL_GENE_NAME].isin(genes_u))
        all_data = df[mask].copy()
    else:
        all_data = df.copy()

    n_rows = len(all_data)
    if n_rows == 0:
        return all_data

    trans_table = str.maketrans("tT", "uU")
    sense_cache = {seq: get_antisense_rna(seq).encode() for seq in all_data[ASO_SEQUENCE].unique()}
    all_data["__temp_sense"] = all_data[ASO_SEQUENCE].map(sense_cache)
    all_data["__temp_idx"] = np.arange(n_rows)
    seq_lengths = all_data[ASO_SEQUENCE].str.len().values
    pre_mrna_cache = {
        g: gene_to_data[g].full_mrna.upper().translate(trans_table).encode() for g in genes_u if g in gene_to_data
    }
    _validate_gene_lengths(pre_mrna_cache.keys(), gene_to_data)
    _validate_codons(pre_mrna_cache.keys(), gene_to_data)

    out = _init_outputs(n_rows)

    for gene_name, group in all_data.groupby(CANONICAL_GENE_NAME, observed=True):
        if gene_name not in pre_mrna_cache:
            continue
        locus_info = gene_to_data[gene_name]
        pre_mrna = pre_mrna_cache[gene_name]

        # setup: locate each ASO's target in the pre-mRNA and map the match to a genomic coordinate.
        group_rows = group["__temp_idx"].values
        match_pos = _match_positions(pre_mrna, group["__temp_sense"].values)
        found = match_pos != -1
        if not found.any():
            continue
        hit_rows = group_rows[found]
        hit_pos = match_pos[found]
        genetic_coordinates = _genetic_coordinates(hit_pos, seq_lengths[hit_rows], locus_info)

        assign_target_position(out, hit_rows, hit_pos, len(pre_mrna), locus_info)
        assign_canonical_splice_junction_distances(out, hit_rows, genetic_coordinates, locus_info)
        assign_closest_splice_junction_distance(out, hit_rows, genetic_coordinates, locus_info)
        assign_host_lengths(out, hit_rows, genetic_coordinates, locus_info)
        assign_distance_to_special_codons(out, hit_rows, genetic_coordinates, locus_info)
        assign_mrna_stop_distance(out, hit_rows, genetic_coordinates, locus_info)
        assign_region(out, hit_rows, genetic_coordinates, locus_info)

    all_data[STRUCTURE_SENSE_START] = out.start
    all_data[STRUCTURE_SENSE_START_FROM_END] = out.start_end
    all_data[STRUCTURE_SENSE_LENGTH] = seq_lengths
    all_data[STRUCT_SENSE_IN_EXON] = out.exon
    all_data[STRUCT_SENSE_IN_EXON_NON_EXCLUSIVE] = out.exon_non_exclusive
    all_data[STRUCT_SENSE_IN_INTRON] = out.intron
    all_data[STRUCT_SENSE_IN_UTR] = out.utr3 | out.utr5
    all_data[STRUCT_SENSE_IN_3UTR] = out.utr3
    all_data[STRUCT_SENSE_IN_5UTR] = out.utr5
    all_data[STRUCT_SENSE_IN_CDS] = out.cds
    all_data[STRUCT_SENSE_IN_CDS_NON_EXCLUSIVE] = out.cds_non_exclusive
    all_data[STRUCTURE_SENSE_TYPE] = out.type
    all_data[STRUCTURE_SENSE_START_NORM] = out.start_norm
    all_data[STRUCTURE_SENSE_START_FROM_END_NORM] = out.start_end_norm
    all_data[STRUCTURE_SENSE_DIST_TO_CANONICAL_STOP] = out.dist_canonical_stop
    all_data[STRUCTURE_SENSE_DIST_TO_CLOSEST_STOP] = out.dist_closest_stop
    all_data[STRUCTURE_SENSE_DIST_TO_CANONICAL_START] = out.dist_canonical_start
    all_data[STRUCTURE_SENSE_DIST_TO_CLOSEST_START] = out.dist_closest_start
    all_data[STRUCTURE_SENSE_MRNA_DIST_TO_CANONICAL_STOP] = out.mrna_dist_canonical_stop
    all_data[STRUCTURE_SENSE_MRNA_DIST_TO_CLOSEST_STOP] = out.mrna_dist_closest_stop
    all_data[STRUCTURE_SENSE_DIST_TO_SPLICE_JUNCTION_EXONIC] = out.dist_sj_exonic
    all_data[STRUCTURE_SENSE_DIST_TO_SPLICE_JUNCTION_INTRONIC] = out.dist_sj_intronic
    all_data[STRUCTURE_SENSE_DIST_TO_CLOSEST_SPLICE_JUNCTION] = out.dist_closest_sj
    all_data[STRUCTURE_SENSE_HOST_EXON_LOG_LENGTH] = out.host_exon_log_len
    all_data[STRUCTURE_SENSE_HOST_INTRON_LOG_LENGTH] = out.host_intron_log_len
    all_data[STRUCTURE_SENSE_JUNCTION_LOGDIST_EXONIC] = out.junction_logdist_exonic
    all_data[STRUCTURE_SENSE_JUNCTION_LOGDIST_INTRONIC] = out.junction_logdist_intronic
    all_data[STRUCTURE_SENSE_JUNCTION_LOGDIST_CLOSEST] = out.junction_logdist_closest

    all_data.drop(columns=["__temp_idx", "__temp_sense"], inplace=True)
    return all_data
