import numpy as np
import pandas as pd

from notebooks.data.OligoAI.parse_chemistry import populate_chemistry
from notebooks.notebook_utils import log_correction, read_cached_gene_to_data
from tauso.data.consts import *
from tauso.util import get_antisense


def get_unique_genes(df):
    """
    Extracts unique genes from a dataframe, excluding blacklist items.
    """
    genes = df[CANONICAL_GENE].unique()
    genes_to_remove = {"HBV", "negative_control", np.nan}

    # Filter using list comprehension
    genes_clean = [g for g in genes if g not in genes_to_remove]
    return genes_clean


def get_filtered_human_data(df):
    """
    Filters the raw DataFrame for Human entries and removes rows with missing Inhibition data.
    Returns a copy to avoid SettingWithCopy warnings.
    """
    # Filter: Organism is Human AND Inhibition is not NaN
    condition = (df[CELL_LINE_ORGANISM] == "human") & (df[INHIBITION].notna())
    return df[condition].copy()


def process_row(row, gene_to_data):
    gene_name = row[CANONICAL_GENE]
    if gene_name not in gene_to_data:
        return -1, 0, ""

    locus_info = gene_to_data[gene_name]
    pre_mrna = locus_info.full_mrna
    antisense = getattr(row, SEQUENCE)
    sense = get_antisense(antisense)
    idx = pre_mrna.find(sense)  # the index in the pre_mrna the sense is found

    return idx, len(antisense), sense


def preprocess_aso_data(csv_path, include_smiles: bool = False):
    """
    Loads raw ASO data, cleans it, calculates log inhibition, and maps sequences.
    """
    if not csv_path.exists():
        raise FileNotFoundError(f"Data file not found at: {csv_path}")

    # 1. Load Data
    all_data = pd.read_csv(str(csv_path), low_memory=False)

    # 2. Filter & Log Correct (Abstracted filtering logic)
    df = get_filtered_human_data(all_data)
    log_correction(df)

    df = df[df[CHEMICAL_PATTERN] != "(S)-cEt/5-methylcytosines/deoxy"].copy()
    df = df[df[CHEMICAL_PATTERN] != "LNA/deoxy"].copy()

    # 3. Filter Specific Genes (using the clean DF)
    genes_clean = get_unique_genes(df)
    gene_to_data = read_cached_gene_to_data(genes_clean)
    df = df[df[CANONICAL_GENE].isin(genes_clean)].copy()

    # 4. Optional: Drop SMILES
    if not include_smiles and SMILES in df.columns:
        df.drop(columns=[SMILES], inplace=True)

    # 5. Sequence Mapping
    # Initialize columns
    df[SENSE_SEQUENCE] = ""
    df[SENSE_START] = -1
    df[SENSE_LENGTH] = 0

    results = df.apply(lambda row: process_row(row, gene_to_data), axis=1)
    df[[SENSE_START, SENSE_LENGTH, SENSE_SEQUENCE]] = pd.DataFrame(results.tolist(), index=df.index)

    # 6. Final Cleanup
    valid_data = df[df[SENSE_START] != -1].copy()

    # 7. Add standard cell line column
    if not CELL_LINE_DEPMAP_PROXY in valid_data.columns:
        print("Adding cell line depmap name proxy")
        valid_data[CELL_LINE_DEPMAP_PROXY] = valid_data[CELL_LINE].map(CELL_LINE_TO_DEPMAP_PROXY_DICT)

    # 8. Add Depmap cell lines
    if not CELL_LINE_DEPMAP in valid_data.columns:
        print("Adding a depmap column")
        valid_data[CELL_LINE_DEPMAP] = valid_data[CELL_LINE_DEPMAP_PROXY].map(CELL_LINE_TO_DEPMAP)

    print(f"Preprocessing complete. Final valid rows: {len(valid_data)}")
    return valid_data


def process_oligo_data(data, min_cohort_size=1, min_cell_line_asos=1, strict_gapmer_patterns=False, verbose=True):
    """
    Takes a raw ASO dataframe, standardizes formats, and sequentially filters out:
    0. Base exclusions (unsupported chemistry, steric blocking, ambiguous genes, missing inhibition, missing cell line)
'm     1. Unmapped sequences
    2. Optional: Non-standard gapmer patterns (filters strictly for 5-10-5 MOE and 3-10-3 cEt)
    3. Sparse cohorts
    4. Sparse cell lines
    """
    if not SENSE_START in data.columns:
        raise ValueError(f"Need {SENSE_START} in data to filter properly! Add the features")

    # Create a working copy to avoid mutating the original dataframe
    data = data.copy()

    # Apply chemistry parsing
    data = populate_chemistry(data)

    initial_total_rows = len(data)

    # ---------------------------------------------------------
    # FILTER: Chemistry (Whitelist MOE and cEt)
    # ---------------------------------------------------------
    rows_before_chem = len(data)
    valid_chemistries = ["MOE/5-methylcytosines/deoxy", "cEt/5-methylcytosines/deoxy"]
    # Assuming MODIFICATION is imported from your consts
    data = data[data[MODIFICATION].isin(valid_chemistries)].copy()
    elim_chemistry = rows_before_chem - len(data)

    # 1. Rename columns to standard consts
    rename_scheme = {
        "aso_sequence_5_to_3": SEQUENCE,
        "Canonical Gene Name": CANONICAL_GENE,
        "cell_line": CELL_LINE,
        "cell_line_species": CELL_LINE_ORGANISM,
        "inhibition_percent": INHIBITION,
        "dosage": VOLUME,
    }
    data = data.rename(columns=rename_scheme)

    # 2. Standardize Cell Line Naming
    cell_line_fixes = {
        "A-431": "A431",
        "A-549": "A549",
        "A459": "A549",  # Typo
        "HepB3": "Hep3B",  # Typo
        "T-24": "T24",
        "SH-SY-5Y": "SH-SY5Y",
        "HuVEC": "HUVEC",  # Capitalization standard
        "hSKM": "hSKMc",  # Consolidate skeletal muscle variants
        "iCell cardiomyocytes (R1017)": "iCell cardiomyocytes2",
    }
    data[CELL_LINE] = data[CELL_LINE].replace(cell_line_fixes)

    # ---------------------------------------------------------
    # FILTER 0: Base Filtering
    # ---------------------------------------------------------
    rows_before = len(data)
    data = data[data["steric_blocking"] == False].copy()
    elim_steric = rows_before - len(data)

    rows_before = len(data)
    data = data[~data[CANONICAL_GENE].str.contains(";", na=False)].copy()
    elim_multigene = rows_before - len(data)

    rows_before = len(data)
    data = data[data[INHIBITION].notna()].copy()
    elim_missing_inhib = rows_before - len(data)

    # Explicitly drop and track missing cell lines
    rows_before = len(data)
    data = data[data[CELL_LINE].notna()].copy()
    elim_missing_cell_line = rows_before - len(data)

    # ---------------------------------------------------------
    # FILTER 1: Unmapped Sequences (SENSE_START == -1)
    # ---------------------------------------------------------
    initial_rows_unmapped = len(data)
    unmapped_series = data[data[SENSE_START] == -1][CANONICAL_GENE].value_counts()
    data = data[data[SENSE_START] != -1].copy()
    elim_unmapped_rows = initial_rows_unmapped - len(data)

    # ---------------------------------------------------------
    # FILTER 2: Strict Gapmer Patterns (OPTIONAL)
    # ---------------------------------------------------------
    elim_patterns = 0
    if strict_gapmer_patterns:
        rows_before_patterns = len(data)

        # 5-10-5 MOE are 20-mers | 3-10-3 cEt are 16-mers
        is_5105_moe = data[CHEMICAL_PATTERN] == "MMMMMddddddddddMMMMM"
        is_3103_cet = data[CHEMICAL_PATTERN] == "CCCddddddddddCCC"

        data = data[is_5105_moe | is_3103_cet].copy()
        elim_patterns = rows_before_patterns - len(data)

    # ---------------------------------------------------------
    # FILTER 3: Sparse Cohorts (< min_cohort_size)
    # ---------------------------------------------------------
    data["cohort_id"] = data[CANONICAL_GENE].astype(str).str.strip() + "_" + data[CELL_LINE].astype(str).str.strip()
    cohort_counts = data["cohort_id"].value_counts()

    valid_cohorts = cohort_counts[cohort_counts >= min_cohort_size].index
    dropped_cohort_series = cohort_counts[cohort_counts < min_cohort_size]

    initial_rows_cohort = len(data)
    initial_cohorts = len(cohort_counts)

    data = data[data["cohort_id"].isin(valid_cohorts)].copy()

    elim_cohort_rows = initial_rows_cohort - len(data)
    final_cohorts = data["cohort_id"].nunique()

    # ---------------------------------------------------------
    # FILTER 4: Sparse Cell Lines (< min_cell_line_asos)
    # ---------------------------------------------------------
    cell_counts = data[CELL_LINE].value_counts()
    valid_cell_lines = cell_counts[cell_counts >= min_cell_line_asos].index
    dropped_cell_series = cell_counts[cell_counts < min_cell_line_asos]

    initial_rows_cell_line = len(data)

    data = data[data[CELL_LINE].isin(valid_cell_lines)].copy()

    elim_cell_rows = initial_rows_cell_line - len(data)

    # ---------------------------------------------------------
    # 5. Map Cell Lines to DepMap Contexts
    # ---------------------------------------------------------
    data[CELL_LINE_DEPMAP_PROXY] = data[CELL_LINE].map(CELL_LINE_TO_DEPMAP_PROXY_DICT)
    data[CELL_LINE_DEPMAP] = data[CELL_LINE_DEPMAP_PROXY].map(CELL_LINE_TO_DEPMAP)

    # ---------------------------------------------------------
    # 6. Print Summary Report
    # ---------------------------------------------------------
    if verbose:
        print("-" * 60)
        print("PROCESSING FILTERING REPORT")
        print("-" * 60)
        print(f"Initial raw rows loaded: {initial_total_rows:,}\n")

        print(f"[0. BASE FILTERING]")
        print(f"Unsupported chemistry (Mixmers/DNA/None): {elim_chemistry:,}")
        print(f"Steric blocking (True) eliminated: {elim_steric:,}")
        print(f"Multiple genes (';' present) eliminated: {elim_multigene:,}")
        print(f"Missing inhibition (NaN) eliminated: {elim_missing_inhib:,}")
        print(f"Missing cell line (NaN) eliminated: {elim_missing_cell_line:,}")

        if SENSE_START in data.columns:
            print(f"\n[1. UNMAPPED SEQUENCES ({SENSE_START} == -1)]")
            print(f"Samples eliminated: {elim_unmapped_rows:,}")

        if strict_gapmer_patterns:
            print(f"\n[2. STRICT GAPMER PATTERNS (5-10-5 MOE / 3-10-3 cEt)]")
            print(f"Non standard sequences eliminated: {elim_patterns:,}")

        print(f"\n[3. COHORT FILTERING (>= {min_cohort_size} samples)]")
        print(f"Cohorts: {initial_cohorts:,} -> {final_cohorts:,} ({len(dropped_cohort_series):,} eliminated)")
        print(
            f"Samples: {initial_rows_cohort:,} -> {initial_rows_cohort - elim_cohort_rows:,} ({elim_cohort_rows:,} eliminated)"
        )

        print(f"\n[4. SPARSE CELL LINE FILTERING (>= {min_cell_line_asos} samples)]")
        print(
            f"Cell Lines: {len(cell_counts):,} -> {len(valid_cell_lines):,} ({len(dropped_cell_series):,} eliminated)"
        )
        print(f"Samples: {initial_rows_cell_line:,} -> {len(data):,} ({elim_cell_rows:,} eliminated)")

        print(f"\nFINAL DATASET: {len(data):,} ASOs")

        print("\n" + "=" * 60)
        print("ELIMINATED GROUPS BREAKDOWN")
        print("=" * 60)

        if SENSE_START in data.columns:
            print(
                f"\n[ELIMINATED UNMAPPED SAMPLES] - {elim_unmapped_rows:,} samples across {len(unmapped_series)} genes:"
            )
            if not unmapped_series.empty:
                for gene, size in unmapped_series.items():
                    print(f"  • {gene}: {size} samples")
            else:
                print("  None")

        print(f"\n[ELIMINATED COHORTS (< {min_cohort_size} samples)] - {len(dropped_cohort_series)} cohorts:")
        if not dropped_cohort_series.empty:
            for cohort, size in list(dropped_cohort_series.items())[:10]:
                print(f"  • {cohort}: {size} samples")
            if len(dropped_cohort_series) > 10:
                print(f"  ... and {len(dropped_cohort_series) - 10} more cohorts.")
        else:
            print("  None")

        print(f"\n[ELIMINATED CELL LINES (< {min_cell_line_asos} samples)] - {len(dropped_cell_series)} cell lines:")
        if not dropped_cell_series.empty:
            for c_line, size in dropped_cell_series.items():
                print(f"  • {c_line}: {size} samples")
        else:
            print("  None")
        print("-" * 60)

    return data
