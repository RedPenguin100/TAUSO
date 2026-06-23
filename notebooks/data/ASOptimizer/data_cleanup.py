"""Clean the raw ASOptimizer dataset into the processed dataset.

Reads ``ORIGINAL_CSV`` (raw ASOptimizer source columns) and writes
``PROCESSED_CSV``. Stages: canonical-gene mapping, cell-line annotation, row
filters (missing inhibition, non-human, non-phosphorothioate, negative
controls), a three-level aggregation (replicates -> cross-probe -> density),
a chemistry whitelist, and dropping modification-scan cohorts.

Usage:
    python -m notebooks.data.ASOptimizer.data_cleanup
"""
import logging

import pandas as pd
from notebooks.data.ASOptimizer.consts import (
    CELL_LINE_ASOPT,
    CHEMICAL_PATTERN_ASOPT,
    DENSITY_ASOPT,
    INHIBITION_ASOPT,
    ISIS_ASOPT,
    LINKAGE_ASOPT,
    LINKAGE_LOCATION_ASOPT,
    MODIFICATION_ASOPT,
    ORIGINAL_CSV,
    PRIMER_PROBE_ASOPT,
    PROCESSED_CSV,
    SEQUENCE_ASOPT,
    TARGET_GENE_ASOPT,
    TRANSFECTION_ASOPT,
    TREATMENT_PERIOD_ASOPT,
    VOLUME_ASOPT,
)

from tauso.data.consts import (
    CANONICAL_GENE_NAME,
    CELL_LINE_DEPMAP,
    CELL_LINE_ORGANISM,
    CELL_LINE_TO_DEPMAP,
    HUMAN,
    cell_line_to_organism,
)

logger = logging.getLogger(__name__)

COHORT = "Cohort"
MOD_SCAN_THRESHOLD = 0.2
VALID_CHEMISTRIES = ["cEt/5-methylcytosines/deoxy", "MOE/5-methylcytosines/deoxy", "LNA/deoxy"]

GENE_TO_CANONICAL = {
    "K-RAS": "KRAS",
    "YAP1": "YAP1",
    "PKK": "KLKB1",
    "SNCA_LNA": "SNCA",
    "ANGPTL2_LNA": "ANGPTL2",
    "TAU": "MAPT",
    "SOD-1": "SOD1",
}

# Aggregation hierarchy: each level drops one grouping column and averages over it.
REPLICATE_GROUP_COLS = [
    ISIS_ASOPT,
    TRANSFECTION_ASOPT,
    LINKAGE_ASOPT,
    LINKAGE_LOCATION_ASOPT,
    DENSITY_ASOPT,
    PRIMER_PROBE_ASOPT,
    CANONICAL_GENE_NAME,
    CELL_LINE_ASOPT,
    VOLUME_ASOPT,
    TREATMENT_PERIOD_ASOPT,
    SEQUENCE_ASOPT,
    MODIFICATION_ASOPT,
    CHEMICAL_PATTERN_ASOPT,
]
CROSS_PROBE_GROUP_COLS = [c for c in REPLICATE_GROUP_COLS if c != PRIMER_PROBE_ASOPT]
FINAL_GROUP_COLS = [c for c in CROSS_PROBE_GROUP_COLS if c != DENSITY_ASOPT]


def load_raw():
    logger.info("loading %s", ORIGINAL_CSV)
    df = pd.read_csv(str(ORIGINAL_CSV), low_memory=False)
    logger.info("loaded %d rows", len(df))
    return df


def assign_canonical_gene(df):
    df = df.copy()
    df[TARGET_GENE_ASOPT] = df[TARGET_GENE_ASOPT].str.upper()
    df[CANONICAL_GENE_NAME] = df[TARGET_GENE_ASOPT].map(GENE_TO_CANONICAL).fillna(df[TARGET_GENE_ASOPT])
    return df


def annotate_cell_line(df):
    df = df.copy()
    df[CELL_LINE_ORGANISM] = df[CELL_LINE_ASOPT].map(cell_line_to_organism)
    df[CELL_LINE_DEPMAP] = df[CELL_LINE_ASOPT].map(CELL_LINE_TO_DEPMAP)
    df["true_length_of_seq"] = df[SEQUENCE_ASOPT].str.len()
    return df


def filter_missing_inhibition(df):
    before = len(df)
    df = df[~df[INHIBITION_ASOPT].isna()]
    logger.info("filter missing inhibition: removed %d, remained %d", before - len(df), len(df))
    return df


def filter_human(df):
    unknown = df[df[CELL_LINE_ORGANISM].isna()][CELL_LINE_ASOPT].unique()
    if len(unknown):
        logger.warning("cell lines with unknown organism (dropped): %s", list(unknown))
    before = len(df)
    df = df[df[CELL_LINE_ORGANISM] == HUMAN]
    logger.info("filter non-human: removed %d, remained %d", before - len(df), len(df))
    return df


def filter_phosphorothioate(df):
    before = len(df)
    df = df[df[LINKAGE_ASOPT] == "phosphorothioate"]
    logger.info("filter phosphorothioate linkage: removed %d, remained %d", before - len(df), len(df))
    return df


def drop_negative_control(df):
    before = len(df)
    df = df[df[TARGET_GENE_ASOPT] != "NEGATIVE_CONTROL"]
    logger.info("drop negative controls: removed %d, remained %d", before - len(df), len(df))
    return df


def aggregate_replicates(df):
    """Average replicate measurements of the exact same experiment."""
    before = len(df)
    agg = df.groupby(REPLICATE_GROUP_COLS, dropna=False).agg({INHIBITION_ASOPT: ["mean", "std", "count"]}).reset_index()
    agg.columns = REPLICATE_GROUP_COLS + ["mean_inhibition", "std_inhibition", "replicate_count"]
    agg["std_inhibition"] = agg["std_inhibition"].fillna(0)
    logger.info("aggregate replicates: %d -> %d", before, len(agg))
    return agg


def aggregate_cross_probe(df):
    """Average results that differ only by the primer/probe set used."""
    before = len(df)
    agg = (
        df.groupby(CROSS_PROBE_GROUP_COLS, dropna=False)
        .agg({"mean_inhibition": ["mean", "std", "count"], "replicate_count": "sum"})
        .reset_index()
    )
    agg.columns = CROSS_PROBE_GROUP_COLS + ["mean_inhibition", "probe_std", "probe_count", "total_replicate_count"]
    agg["probe_std"] = agg["probe_std"].fillna(0)
    logger.info("aggregate cross-probe: %d -> %d", before, len(agg))
    return agg


def aggregate_density(df):
    """Average results that differ only by cell density; keep densities as a list."""
    before = len(df)
    agg = (
        df.groupby(FINAL_GROUP_COLS, dropna=False)
        .agg(
            {
                "mean_inhibition": "mean",
                "probe_std": "mean",
                "probe_count": "sum",
                "total_replicate_count": "sum",
                DENSITY_ASOPT: lambda x: ", ".join(x.astype(str).unique()),
            }
        )
        .reset_index()
    )
    agg = agg.rename(columns={"mean_inhibition": INHIBITION_ASOPT})
    logger.info("aggregate density: %d -> %d", before, len(agg))
    return agg


def filter_chemistry(df):
    """Keep only the whitelisted pure chemistries; drop mixmers and sparse classes."""
    before = len(df)
    df = df[df[MODIFICATION_ASOPT].isin(VALID_CHEMISTRIES)]
    logger.info("filter chemistry: removed %d, remained %d", before - len(df), len(df))
    return df


def add_cohort(df):
    df = df.copy()
    df[COHORT] = df[CANONICAL_GENE_NAME] + " (" + df[CELL_LINE_ASOPT] + ")"
    return df


def compute_cohort_diversity(df):
    """Per-cohort sequence-diversity stats (unique sequences / unique molecular designs).

    Requires the COHORT column to be present (see add_cohort).
    """
    unique_designs = df.groupby(
        [COHORT, SEQUENCE_ASOPT, MODIFICATION_ASOPT, CHEMICAL_PATTERN_ASOPT, LINKAGE_ASOPT, LINKAGE_LOCATION_ASOPT]
    ).size().reset_index()
    cohort_designs = unique_designs.groupby(COHORT).size().to_frame("Unique_Molecular_Designs")
    cohort_seqs = df.groupby(COHORT)[SEQUENCE_ASOPT].nunique().to_frame("Unique_Sequences")
    diversity = cohort_seqs.join(cohort_designs).reset_index()
    diversity["Diversity_Ratio"] = diversity["Unique_Sequences"] / diversity["Unique_Molecular_Designs"]
    return diversity


def drop_mod_scan_cohorts(df, threshold=MOD_SCAN_THRESHOLD):
    """Drop modification-scan cohorts: those whose sequence-diversity ratio is below threshold."""
    df = add_cohort(df)
    diversity = compute_cohort_diversity(df)
    to_remove = diversity[diversity["Diversity_Ratio"] < threshold][COHORT].tolist()
    before = len(df)
    df = df[~df[COHORT].isin(to_remove)]
    logger.info("modification-scan cohorts dropped (%d): %s", len(to_remove), to_remove)
    logger.info("drop mod-scan cohorts: removed %d, remained %d", before - len(df), len(df))
    return df


def finalize(df):
    df = df.copy()
    df = df.rename(columns={"mean_inhibition": INHIBITION_ASOPT})
    df["index_v2"] = range(len(df))
    df[CELL_LINE_ORGANISM] = HUMAN
    return df


def run_pipeline(through_mod_scan: bool = True):
    """Run the full cleaning pipeline. If through_mod_scan is False, stop after the
    cross-probe aggregation -- the stage the diversity plot visualizes."""
    df = load_raw()
    df = assign_canonical_gene(df)
    df = annotate_cell_line(df)
    df = filter_missing_inhibition(df)
    df = filter_human(df)
    df = filter_phosphorothioate(df)
    df = drop_negative_control(df)
    df = aggregate_replicates(df)
    df = aggregate_cross_probe(df)
    if not through_mod_scan:
        return df
    df = drop_mod_scan_cohorts(df)
    df = finalize(df)
    return df


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s | %(message)s")
    df = run_pipeline()
    logger.info("saving %d rows -> %s", len(df), PROCESSED_CSV)
    df.to_csv(PROCESSED_CSV)


if __name__ == "__main__":
    main()
