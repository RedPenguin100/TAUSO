import logging

from notebooks.data.OligoAI.parse_chemistry import assign_chemistry
from tauso.data.consts import *
from tauso.util import get_antisense

logger = logging.getLogger(__name__)


def process_row(row, gene_to_data):
    gene_name = row[CANONICAL_GENE_NAME]
    if gene_name not in gene_to_data:
        return -1, 0, ""

    locus_info = gene_to_data[gene_name]
    pre_mrna = locus_info.full_mrna
    antisense = getattr(row, ASO_SEQUENCE)
    sense = get_antisense(antisense)
    idx = pre_mrna.find(sense)  # the index in the pre_mrna the sense is found

    return idx, len(antisense), sense


def process_oligo_data_rename(data):
    """Translate the raw OligoAI source columns to TAUSO's canonical names.
    Identity entries (cell_line, inhibition_percent) are intentionally omitted -- the
    raw and canonical names already match after the consts.py standardization, so the
    rename map only carries the columns that actually need renaming.
    """
    rename_scheme = {
        "aso_sequence_5_to_3": ASO_SEQUENCE,
        "Canonical Gene Name": CANONICAL_GENE_NAME,
        "cell_line_species": CELL_LINE_ORGANISM,
        "dosage": VOLUME_NM,
        "cells_per_well": DENSITY_CELLS_PER_WELL,
        "transfection_method": TRANSFECTION_RAW,
    }
    data = data.rename(columns=rename_scheme)
    return data


CELL_LINE_FIXES = {
    "A-431": "A431",
    "A-549": "A549",
    "A459": "A549",  # typo
    "HepB3": "Hep3B",  # typo
    "T-24": "T24",
    "SH-SY-5Y": "SH-SY5Y",
    "HuVEC": "HUVEC",
    "hSKM": "hSKMc",  # consolidate skeletal muscle variants
    "iCell cardiomyocytes (R1017)": "iCell cardiomyocytes2",
}
# Single high-affinity sugar (MOE / cEt / LNA) or plain DNA. LNA is supported even though the
# ASO Atlas contains none. Mixed-sugar oligos ("mixmer", >1 high-affinity sugar) are excluded.
SUPPORTED_CHEMISTRIES = [
    "MOE/5-methylcytosines/deoxy",
    "cEt/5-methylcytosines/deoxy",
    "LNA/5-methylcytosines/deoxy",
    "DNA",
]
STRICT_GAPMER_PATTERNS = ["MMMMMddddddddddMMMMM", "CCCddddddddddCCC"]  # 5-10-5 MOE, 3-10-3 cEt


def _keep(data, mask, reason):
    """Keep rows where mask is True, logging how many were dropped (including 0, so every filter is visible)."""
    kept = data[mask].copy()
    logger.info("  dropped %7d  %s", len(data) - len(kept), reason)
    return kept


def _standardize(data):
    data = assign_chemistry(data)
    data = process_oligo_data_rename(data)
    data[CELL_LINE] = data[CELL_LINE].replace(CELL_LINE_FIXES)
    return data


def _filter_supported_chemistry(data):
    return _keep(data, data[MODIFICATION_STRING].isin(SUPPORTED_CHEMISTRIES),
                 "mixed-sugar chemistry (mixmer)")


def _filter_valid_target(data):
    data = _keep(data, data["steric_blocking"] == False, "steric blocking")
    data = _keep(data, data[CANONICAL_GENE_NAME].notna(), "missing canonical gene")
    data = _keep(data, ~data[CANONICAL_GENE_NAME].str.contains(";", na=False), "multi-gene target")
    data = _keep(data, data[INHIBITION_PERCENT].notna(), "missing inhibition")
    data = _keep(data, data[CELL_LINE].notna(), "missing cell line")
    return data


def _filter_mapped(data):
    return _keep(data, data[STRUCTURE_SENSE_START] != -1, "sequence unmapped to target")


def _filter_strict_gapmers(data):
    return _keep(data, data[CHEMICAL_PATTERN].isin(STRICT_GAPMER_PATTERNS), "non-strict gapmer pattern")


def _filter_sparse(data, min_cohort_size, min_cell_line_asos):
    data["cohort_id"] = (data[CANONICAL_GENE_NAME].astype(str).str.strip() + "_"
                         + data[CELL_LINE].astype(str).str.strip())
    cohort_counts = data["cohort_id"].value_counts()
    data = _keep(data, data["cohort_id"].isin(cohort_counts[cohort_counts >= min_cohort_size].index),
                 f"cohort with < {min_cohort_size} ASOs")
    cell_counts = data[CELL_LINE].value_counts()
    data = _keep(data, data[CELL_LINE].isin(cell_counts[cell_counts >= min_cell_line_asos].index),
                 f"cell line with < {min_cell_line_asos} ASOs")
    return data


def process_oligo_data(data, min_cohort_size=1, min_cell_line_asos=1, strict_gapmer_patterns=False):
    """Standardize and filter the raw ASO table down to the analysis-ready set.

    Requires STRUCTURE_SENSE_START (run structure features first). Each filter logs its drop count.
    """
    if STRUCTURE_SENSE_START not in data.columns:
        raise ValueError(f"Need {STRUCTURE_SENSE_START} to filter; compute structure features first.")

    data = _standardize(data.copy())
    logger.info("process_oligo_data: %d rows in", len(data))

    data = _filter_supported_chemistry(data)
    data = _filter_valid_target(data)
    data = _filter_mapped(data)
    if strict_gapmer_patterns:
        data = _filter_strict_gapmers(data)
    data = _filter_sparse(data, min_cohort_size, min_cell_line_asos)

    data[CELL_LINE_DEPMAP_PROXY] = data[CELL_LINE].map(resolve_depmap_proxy)
    data[CELL_LINE_DEPMAP] = data[CELL_LINE].map(resolve_depmap_id)
    logger.info("process_oligo_data: %d rows out", len(data))
    return data
