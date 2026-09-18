import logging
import os
from collections import namedtuple

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

logger = logging.getLogger(__name__)

from ...data.data import get_data_dir
from ...frames import release_arrow_memory

# Cell-line tokens that carry no usable identity. They skip the exact (Tier-1)
# lookup and resolve to the gene-level estimate, so a missing/unknown cell line
# never matches TTDB's catch-all "Unknown" cell type by accident.
MISSING_CELL_TOKENS = {"", "nan", "none", "na", "unknown", "generic"}

# Result of a half-life query. `source` records which tier supplied the value
# (cell-specific measurement, gene-level geometric mean, or absent -> NaN).
HalfLifeResult = namedtuple("HalfLifeResult", ["half_life", "source"])

# Columns the loader needs from the raw TTDB dataset. The full file (21 columns,
# distributed as csv.gz on Zenodo) is converted to a lean Parquet of just these
# columns by `tauso setup-mrna-halflife`.
#
# `gene_name_x` / `gene_name_y` are pandas merge() suffixes present in TTDB's
# officially distributed file (species_stability_no_threshold.csv.gz, the Google
# Drive mirror linked from https://sysbio.gzzoc.com/ttdb/download.html) — TTDB's
# standardization pipeline joined two tables that each carry a `gene_name`
# column. They are NOT created by this repo. Verified across all 2.17M rows:
# the SAME gene per row (identical after upper()+strip()); the only raw
# difference (~35% of rows) is letter case, with `gene_name_y` the uppercase-
# normalized form (~99.8% all-upper vs ~64% for `_x`). We use `_y`.
HALFLIFE_SOURCE_COLUMNS = [
    "gene_name_y",
    "species_name",
    "cell_type",
    "half_life",
    "condition",
    "r_squared",
]

# TTDB is multi-species; the dataset is human throughout.
HALFLIFE_SPECIES = "Human"

# A study's baseline arm carries whatever name that study gave it; the rest are treatments.
HALFLIFE_BASELINE_CONDITIONS = (
    "WT",
    "uninfected",
    "DRB release control",
    "FP control",
    "ZAK control sham vs UVB",
)


def _distinct(series):
    """The column's distinct values: its categories when it has them, its values otherwise.

    A categorical column already holds each distinct value once, and the half-life table has
    a few dozen of them over millions of rows. Matching against those rather than against a
    string per row is what keeps the filter from expanding the column it is reading.
    """
    values = series.cat.categories if isinstance(series.dtype, pd.CategoricalDtype) else series.unique()
    return pd.Index(values)


def _select(df, column, wanted, absent_message):
    """Rows whose `column`, stripped, is one of `wanted`, ignoring case."""
    distinct = _distinct(df[column])
    stripped = distinct.astype(str).str.strip()
    wanted = {str(w).strip().casefold() for w in wanted}
    hit = df[column].isin(distinct[stripped.str.casefold().isin(wanted)])
    if not hit.any():
        raise ValueError(f"{absent_message} in the half-life data; it has {sorted(stripped.unique())}.")
    return df[hit]


def select_species(df, species):
    """Rows for `species`. Raises if absent, naming what the file has."""
    return _select(df, "species_name", [species], f"No {species!r} rows")


def select_conditions(df, conditions):
    """Rows for the baseline `conditions`. Raises if none match, naming what the file has."""
    return _select(df, "condition", conditions, f"No {sorted(conditions)} rows")


def load_halflife_mapping(species=HALFLIFE_SPECIES, conditions=HALFLIFE_BASELINE_CONDITIONS):
    """
    Loads the TTDB data using the tauso directory structure.
    Returns a dictionary: {(Gene_Symbol, Cell_Line): Half_Life_Hours}
    """
    data_dir = get_data_dir()
    path = os.path.join(data_dir, "mrna_half_life.parquet")

    logger.info(f"Loading half-life data from {path}...")

    if not os.path.exists(path):
        raise FileNotFoundError(f"Data file not found at {path}. Run 'tauso setup-mrna-halflife' first.")

    # 1. Load necessary columns. Parquet is columnar, so this reads only these
    # (we use 'gene_name_y' standardized and 'condition' for filtering).
    #
    # The text columns are read as categories: over 2.17M rows they hold 5 species, 26
    # conditions, 46 cell types and 63k genes, so a category code per row and one copy of
    # each distinct string is a tenth of the size of a Python string per row.
    have = set(pq.ParquetFile(path).schema_arrow.names)
    missing = [c for c in HALFLIFE_SOURCE_COLUMNS if c not in have]
    if missing:
        raise ValueError(f"{path} is missing {missing}. Rebuild it with 'tauso setup-mrna-halflife --force'.")
    df = pq.read_table(path, columns=HALFLIFE_SOURCE_COLUMNS).to_pandas(strings_to_categorical=True)

    # 2. Rename for clarity
    df = df.rename(columns={"gene_name_y": "gene", "cell_type": "cell_line"})

    # 3. One species only.
    df = select_species(df, species)

    # 4. Baseline arms only, not stress responses.
    df = select_conditions(df, conditions)

    # 5. Quality Control Filter (Smart R_Squared)
    # Convert to numeric, forcing errors to NaN
    df["r_squared"] = pd.to_numeric(df["r_squared"], errors="coerce")

    # LOGIC:
    # - If R^2 is >= 0.7: KEEP (Good fit)
    # - If R^2 is NaN: KEEP (Benefit of the doubt / source didn't report it)
    # - If R^2 is < 0.7: DROP (Proven bad fit / failed experiment)
    df = df[(df["r_squared"] >= 0.7) | (df["r_squared"].isna())]

    # 6. Numerical Cleaning
    df["half_life"] = pd.to_numeric(df["half_life"], errors="coerce")
    df = df.dropna(subset=["gene", "half_life"])

    # Filter out negative or zero half-lives (physically impossible)
    df = df[df["half_life"] > 0]

    # 7. Clip Artifacts
    # We clip to 48h, but thanks to step 5, we aren't clipping "failed" 1M hour experiments
    df["half_life"] = df["half_life"].clip(upper=48.0)

    # 8. String Normalization
    df["gene"] = df["gene"].str.upper().str.strip()
    df["cell_line"] = df["cell_line"].str.strip()

    # 9. Handle Duplicates using GEOMETRIC MEAN
    # Since we have filtered out the "noisy" experiments, the remaining
    # duplicates are likely valid replicates. Geometric mean is best for rates.
    # gmean == exp(mean(log(x))); compute it vectorized via the native groupby
    # mean. groupby(...).apply(gmean) calls scipy once per group (~330k groups),
    # which is ~450x slower for an identical result.
    df["_log_half_life"] = np.log(df["half_life"])
    df_clean = np.exp(df.groupby(["gene", "cell_line"], observed=True)["_log_half_life"].mean())

    # 10. Convert to Dictionary
    mapping = df_clean.to_dict()

    logger.info(f"Successfully loaded {len(mapping)} specific (Gene+Cell) {species} stability profiles.")
    release_arrow_memory()
    return mapping


class HalfLifeProvider:
    def __init__(self, mapping):
        self.exact_map = mapping

        # Convert mapping (dict) back to DataFrame for easier tier-level aggregation
        # Mapping structure: {(Gene, Cell): Half_Life}
        df_map = pd.DataFrame([{"gene": k[0], "cell": k[1], "hl": v} for k, v in mapping.items()])

        # --- Tier 2 (Gene Stats) ---
        # Gene-level estimate: geometric mean of the gene's clipped half-life
        # across all available cell lines. Computed vectorized as exp(mean(log))
        # to avoid the per-group scipy gmean call.
        df_map["_log_hl"] = np.log(df_map["hl"])
        gene_geom = np.exp(df_map.groupby("gene", observed=True)["_log_hl"].mean())
        self.gene_geom_mean = gene_geom.to_dict()  # gene -> geometric-mean half-life

        logger.info(f"Provider ready. {len(self.gene_geom_mean)} genes with a gene-level estimate.")

    def get_halflife(self, gene, cell_line):
        """
        Resolves a half-life with a two-tier fallback.

        Tier 1: cell-specific TTDB measurement.
        Tier 2: gene-level geometric mean across cell lines.
        Absent: gene not in TTDB -> NaN, so the model treats it as missing
        rather than being handed a near-constant global mean.
        """
        g = str(gene).upper().strip()
        c = str(cell_line).strip()

        # Tier 1: cell-specific measurement (skipped for missing/placeholder cells)
        if c.lower() not in MISSING_CELL_TOKENS and (g, c) in self.exact_map:
            return HalfLifeResult(self.exact_map[(g, c)], "Experimental (Specific)")

        # Tier 2: gene-level geometric mean across cell lines
        if g in self.gene_geom_mean:
            return HalfLifeResult(self.gene_geom_mean[g], "Imputed (Gene GeomMean)")

        # Absent: gene unknown to TTDB
        return HalfLifeResult(np.nan, "Absent")
