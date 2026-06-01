"""Ribo-seq features per ASO target site.

Two GWIPS-viz tracks from Wagner et al. 2020 (HEK293T, Sel-TCP-seq, GSE139131):
40S unselected scanning footprints and 80S unselected elongation footprints.
Same cell line for both. See ``tauso setup-riboseq``.

Column prefixes: ``ribo_40s_*`` / ``ribo_80s_*``.
"""

import logging
import time
from pathlib import Path

import numpy as np
import pyBigWig

from ...data.consts import CANONICAL_GENE, STRUCTURE_SENSE_LENGTH, STRUCTURE_SENSE_START
from ...data.data import get_data_dir

logger = logging.getLogger(__name__)


_TRACKS = {
    "40s": ("human_unselected_40S.RiboProElong.bw", "ribo_40s"),  # scanning / pre-initiation
    "80s": ("human_unselected_80S.RiboCov.bw", "ribo_80s"),  # elongation footprints
}


def get_ribo_bigwig_path(track="40s"):
    """Resolve the path to a ribo-seq bigWig in TAUSO_DATA_DIR.

    track="40s" returns the Wagner20 40S scanning track (back-compat default).
    track="80s" returns the Wagner20 80S elongation footprint track.
    """
    if track not in _TRACKS:
        raise ValueError(f"unknown ribo-seq track {track!r}; choose from {sorted(_TRACKS)}")
    filename, _prefix = _TRACKS[track]
    bw_path = Path(get_data_dir()) / filename
    if not bw_path.exists():
        raise FileNotFoundError(f"Ribo-seq bigWig not found at {bw_path}. Run 'tauso setup-riboseq' to download it.")
    return str(bw_path)


def get_ribo_40s_human_data():
    """Deprecated: use ``get_ribo_bigwig_path("40s")``."""
    return get_ribo_bigwig_path("40s")


def get_feature_prefix(track):
    """Return the column-name prefix for a given track."""
    return _TRACKS[track][1]


def feature_names(flanks, how, prefix="ribo_40s"):
    """Return the ordered list of ribo-seq feature column names for given flanks and reduction."""
    names = [f"{prefix}_gene_{how}"]
    for f in flanks:
        names.append(f"{prefix}_f{f}_{how}")
        if f > 0:
            names.append(f"{prefix}_upstream_{f}_{how}")
            names.append(f"{prefix}_downstream_{f}_{how}")
    return names


def reduce_values(values, how):
    if values.size == 0:
        return 0.0
    if how == "sum":
        return values.sum()
    elif how == "mean":
        return values.mean()
    elif how == "max":
        return values.max()
    elif how == "nz_mean":
        nz = values[values > 0]
        return nz.mean() if nz.size > 0 else 0.0
    elif how == "nz_frac":
        return float((values > 0).sum()) / values.size
    raise ValueError(f"No reducing logic for value: {how}")


def _query_bw(bw, chrom, start, end, chrom_len):
    """Query a clamped BigWig region. Returns float array or None for empty/invalid regions."""
    s = max(0, start)
    e = min(end, chrom_len)
    if e <= s:
        return None
    return np.nan_to_num(np.array(bw.values(chrom, s, e), dtype=float), nan=0.0)


def _row_target_features(ts, te, strand, bw, chrom, chrom_len, flanks, how, prefix="ribo_40s"):
    """Compute all flank-based features for a single target position (ts, te)."""
    feat = {}
    for f in flanks:
        arr = _query_bw(bw, chrom, ts - f, te + f, chrom_len)
        feat[f"{prefix}_f{f}_{how}"] = reduce_values(arr, how) if arr is not None else 0.0
        if f > 0:
            arr_left = _query_bw(bw, chrom, ts - f, ts, chrom_len)
            arr_right = _query_bw(bw, chrom, te, te + f, chrom_len)
            val_left = reduce_values(arr_left, how) if arr_left is not None else 0.0
            val_right = reduce_values(arr_right, how) if arr_right is not None else 0.0
            if strand == "+":
                feat[f"{prefix}_upstream_{f}_{how}"] = val_left
                feat[f"{prefix}_downstream_{f}_{how}"] = val_right
            else:
                feat[f"{prefix}_upstream_{f}_{how}"] = val_right
                feat[f"{prefix}_downstream_{f}_{how}"] = val_left
    return feat


def process_gene_group(bw_path, gene_rows, flanks, how, prefix="ribo_40s"):
    """
    Process all rows for one gene using a private BigWig handle (thread-safe).

    The gene-level BigWig query is made once and shared across all rows.
    Returns {original_df_index: feature_dict}.
    Raises KeyError(chrom) when the contig is absent from the BigWig file.
    """
    bw = pyBigWig.open(bw_path)
    try:
        first = gene_rows.iloc[0]
        chrom = first["chrom"]
        chroms = bw.chroms()
        if not isinstance(chrom, str) or chrom not in chroms:
            raise KeyError(chrom)

        chrom_len = chroms[chrom]
        gene_start = int(first["gene_start"])
        gene_end = int(first["gene_end"])

        # Compute gene-level value once for the whole group
        gene_arr = _query_bw(bw, chrom, gene_start, gene_end, chrom_len)
        gene_val = reduce_values(gene_arr, how) if gene_arr is not None else 0.0
        gene_key = f"{prefix}_gene_{how}"

        results = {}
        for idx, row in gene_rows.iterrows():
            ts, te = int(row["target_start"]), int(row["target_end"])
            feat = {gene_key: gene_val}
            feat.update(_row_target_features(ts, te, row["strand"], bw, chrom, chrom_len, flanks, how, prefix))
            results[idx] = feat
        return results
    finally:
        bw.close()


def add_genomic_coordinates(aso_df, mapper):
    """
    Map each ASO row to genomic coordinates.

    Calls the mapper once per unique gene (not once per row), then broadcasts
    the gene-level info and computes target positions with vectorized arithmetic.
    """
    out = aso_df.copy()

    # One mapper call per unique gene
    unique_genes = out[CANONICAL_GENE].dropna().unique()
    t0 = time.perf_counter()
    gene_info = {g: mapper.get_gene_coords(g) for g in unique_genes}
    logger.info(
        "add_genomic_coordinates: %d unique genes looked up in %.2fs",
        len(unique_genes),
        time.perf_counter() - t0,
    )

    def _field(gene, key):
        info = gene_info.get(gene)
        return info[key] if info else np.nan

    out["chrom"] = out[CANONICAL_GENE].map(lambda g: _field(g, "chrom"))
    out["gene_start"] = out[CANONICAL_GENE].map(lambda g: _field(g, "gene_start"))
    out["gene_end"] = out[CANONICAL_GENE].map(lambda g: _field(g, "gene_end"))
    out["strand"] = out[CANONICAL_GENE].map(lambda g: _field(g, "strand"))

    # Vectorized target-position arithmetic (preserves NaN for unmapped genes)
    valid = out["chrom"].notna()
    gene_len = out["gene_end"] - out["gene_start"]
    L = out[STRUCTURE_SENSE_START].astype(float)
    k = out[STRUCTURE_SENSE_LENGTH].astype(float)

    # STRUCTURE_SENSE_START is 0-based (from pre_mrna.find(sense)); gene_start/end are 0-based half-open.
    t_start_plus = out["gene_start"] + L
    t_start_minus = out["gene_start"] + (gene_len - L - k)
    t_start = np.where(out["strand"] == "+", t_start_plus, t_start_minus)

    out["target_start"] = np.where(valid, t_start, np.nan)
    out["target_end"] = np.where(valid, t_start + k, np.nan)

    logger.info(
        "add_genomic_coordinates: %d rows total, %d mapped",
        len(out),
        int(valid.sum()),
    )
    return out
