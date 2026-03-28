import os

import numpy as np
import pandas as pd
from scipy.stats import gmean

from ...data.data import get_data_dir
from ...data.consts import CELL_LINE, CANONICAL_GENE


def load_halflife_mapping():
    """
    Loads the TTDB data using the tauso directory structure.
    Returns a dictionary: {(Gene_Symbol, Cell_Line): Half_Life_Hours}
    """
    data_dir = get_data_dir()
    path = os.path.join(data_dir, "mrna_half_life.csv.gz")

    print(f"Loading half-life data from {path}...")

    if not os.path.exists(path):
        raise FileNotFoundError(f"Data file not found at {path}. Run 'tauso setup-mrna-halflife' first.")

    # 1. Load necessary columns
    # We use 'gene_name_y' (standardized) and 'condition' (for filtering)
    cols_to_use = ["gene_name_y", "cell_type", "half_life", "condition", "r_squared"]

    # Load data
    df = pd.read_csv(path, compression="gzip", usecols=cols_to_use)

    # 2. Rename for clarity
    df = df.rename(columns={"gene_name_y": "gene", "cell_type": "cell_line"})

    # 3. Filter for Wild Type (WT) only
    # Crucial: We only want baseline stability, not stress responses.
    df = df[df["condition"].astype(str).str.strip().str.upper() == "WT"]

    # 4. Quality Control Filter (Smart R_Squared)
    # Convert to numeric, forcing errors to NaN
    df["r_squared"] = pd.to_numeric(df["r_squared"], errors="coerce")

    # LOGIC:
    # - If R^2 is >= 0.7: KEEP (Good fit)
    # - If R^2 is NaN: KEEP (Benefit of the doubt / source didn't report it)
    # - If R^2 is < 0.7: DROP (Proven bad fit / failed experiment)
    df = df[(df["r_squared"] >= 0.7) | (df["r_squared"].isna())]

    # 5. Numerical Cleaning
    df["half_life"] = pd.to_numeric(df["half_life"], errors="coerce")
    df = df.dropna(subset=["gene", "half_life"])

    # Filter out negative or zero half-lives (physically impossible)
    df = df[df["half_life"] > 0]

    # 6. Clip Artifacts
    # We clip to 48h, but thanks to step 4, we aren't clipping "failed" 1M hour experiments
    df["half_life"] = df["half_life"].clip(upper=48.0)

    # 7. String Normalization
    df["gene"] = df["gene"].str.upper().str.strip()
    df["cell_line"] = df["cell_line"].str.strip()

    # 8. Handle Duplicates using GEOMETRIC MEAN
    # Since we have filtered out the "noisy" experiments, the remaining
    # duplicates are likely valid replicates. Geometric mean is best for rates.
    df_clean = df.groupby(["gene", "cell_line"])["half_life"].apply(gmean)

    # 9. Convert to Dictionary
    mapping = df_clean.to_dict()

    print(
        f"Successfully loaded {len(mapping)} specific (Gene+Cell) stability profiles."
    )
    return mapping


class HalfLifeProvider:
    def __init__(self, mapping):
        self.exact_map = mapping

        # Convert mapping (dict) back to DataFrame for easier tier-level aggregation
        # Mapping structure: {(Gene, Cell): Half_Life}
        df_map = pd.DataFrame(
            [{"gene": k[0], "cell": k[1], "hl": v} for k, v in mapping.items()]
        )

        # --- Tier 2 (Gene Stats) ---
        # Paper: "Geometric mean of that clipped gene’s half-life across all available cell-lines"

        # We group by gene and aggregate using gmean (and std/count for metadata)
        # Note: We use ddof=1 for std dev to estimate sample standard deviation
        gene_stats_df = df_map.groupby("gene")["hl"].agg(
            geom_mean=gmean, count="count", std="std"
        )
        # Convert to dict for O(1) lookup: gene -> {stats}
        self.gene_stats = gene_stats_df.to_dict(orient="index")

        # --- Tier 3 (Global Fallback) ---
        # Paper: "Global geometric mean of the clipped dataset"
        all_values = df_map["hl"].values

        if len(all_values) > 0:
            self.global_val = gmean(all_values)
            self.global_std = np.std(all_values, ddof=1)
        else:
            self.global_val = 0.0
            self.global_std = 0.0

        self.global_count = len(all_values)

        print(
            f"Provider ready. Global Geometric Mean: {self.global_val:.2f}h (N={self.global_count})"
        )

    def get_halflife(self, gene, cell_line):
        """
        Retrieves half-life with hierarchical fallback.
        """
        g = str(gene).upper().strip()
        c = str(cell_line).strip()

        # Tier 1: Exact Match (Specific Measurement)
        if (g, c) in self.exact_map:
            return self.exact_map[(g, c)], "Experimental (Specific)", 1, 0.0

        # Tier 2: Gene Estimate (Gene-Level)
        if g in self.gene_stats:
            stats = self.gene_stats[g]
            # RETURNING GEOMETRIC MEAN NOW
            return (
                stats["geom_mean"],
                "Imputed (Gene GeomMean)",
                int(stats["count"]),
                stats["std"],
            )

        # Tier 3: Global Baseline
        return self.global_val, "Imputed (Global GeomMean)", 0, self.global_std


# We map cell lines from the data to their closest relatives in the mRNA half life database.

cell_line_mapping = {
    # --- Direct / Strong Matches ---
    "HepG2": "HepG2",
    "HeLa": "Hela",
    "Hela": "Hela",
    "MCF7": "MCF-7",  # Breast -> Breast
    "K-562": "K562",  # CML -> CML
    "nan": "Unknown",  # Handling missing values
    # --- Liver (Hepatocellular Carcinoma) ---
    # Target: HepG2
    "Hep3B": "HepG2",
    "HepB3": "HepG2",  # Likely typo for Hep3B
    "HepG2/Hep3B": "HepG2",
    "Huh7": "HepG2",
    "HepaRG": "HepG2",
    "SNU-449": "HepG2",
    # --- Prostate ---
    # Target: C4-2 (LNCaP derivative)
    "LNCaP": "C4-2",  # C4-2 is a subline of LNCaP
    "VCaP": "C4-2",  # Prostate cancer proxy
    "PC3": "C4-2",  # Prostate cancer proxy
    # --- Kidney ---
    # Target: HEK293
    "HK-2": "HEK293",  # Kidney Proximal Tubule -> Embryonic Kidney
    # --- Lung ---
    # Target: MRC5VA (Lung Fibroblast)
    # Following your precedent of mapping Lung Cancer -> Lung Fibroblast
    "A549": "MRC5VA",
    "A-549": "MRC5VA",
    "A459": "MRC5VA",  # Likely typo for A549
    "NCI-H460": "MRC5VA",
    # --- Skin / Melanoma ---
    # Target: N/TERT-1_keratinocytes
    "SK-MEL-28": "N/TERT-1_keratinocytes",
    "G-361": "N/TERT-1_keratinocytes",
    "A431": "N/TERT-1_keratinocytes",
    "A-431": "N/TERT-1_keratinocytes",
    # --- Blood / Immune / Leukemia ---
    # Target: K562, U937, T_cells
    "MM.1R": "K562",  # Myeloma -> Leukemia proxy (per your previous logic)
    "THP-1": "U937",  # AML/Monocyte -> Lymphoma/Monocyte (Very similar suspension lines)
    "Jurkat": "T_cells",  # T-cell Leukemia -> T_cells
    "T24": "Hela",  # Bladder Carcinoma -> Hela (Generic Epithelial proxy due to no bladder target)
    "T-24": "Hela",
    # --- Brain / Neuronal ---
    # Target: neural_precursor_cells
    "SH-SY5Y": "neural_precursor_cells",
    "SH-SY-5Y": "neural_precursor_cells",
    "U251": "neural_precursor_cells",
    "A172": "neural_precursor_cells",  # Glioblastoma
    # --- Mesenchymal / Stromal / Muscle / Endothelial ---
    # Target: MRC5VA (Fibroblast as general mesenchymal proxy)
    "hSKMc": "MRC5VA",  # Skeletal Muscle
    "hSKM": "MRC5VA",
    "HuVEC": "MRC5VA",  # Endothelial
    "HUVEC": "MRC5VA",
    "SCA2-04": "MRC5VA",  # Likely patient-derived fibroblast
    "differentiated human adipocytes": "MRC5VA",  # Adipocytes are mesenchymal
    "SW872": "MRC5VA",  # Liposarcoma
    # --- Ovarian ---
    "SKOV3": "Hela",  # Ovarian Adenocarcinoma -> Hela (Generic Epithelial proxy)
    # --- Other / Unmapped High Specificity ---
    # These are highly differentiated and don't map well to cancer lines or fibroblasts.
    "iCell cardiomyocytes2": "Unknown",
    "iCell cardiomyocytes (R1017)": "Unknown",
}


def populate_mrna_halflife_features(all_data, provider):
    """
    Enriches the dataframe with mRNA half-life features by mapping cell lines
    to TTDB proxies and querying the HalfLifeProvider.

    Args:
        all_data (pd.DataFrame): Dataframe containing CANONICAL_GENE and CELL_LINE columns.
                      ['mRNA_HalfLife', 'HalfLife_Source', 'Mapped_Cell_Proxy']
    """
    print(f"Calculating stability features for {len(all_data)} rows...")

    features = ["mRNA_HalfLife", "HalfLife_Source", "Mapped_Cell_Proxy"]

    # 2. Define the logic for a single row
    def _get_halflife_features(row):
        # Extract keys using global constants
        gene = row[CANONICAL_GENE]
        cell = row[CELL_LINE]

        # A. Map to TTDB Proxy
        # Uses the imported 'cell_line_mapping' dict to standardize names
        proxy_cell = cell_line_mapping.get(cell, cell)

        # B. Query the Provider
        # Returns: (half_life, source, n_support, std_dev)
        hl_val, source, n, std = provider.get_halflife(gene, proxy_cell)

        # C. Clip Artifacts (Standardize max duration to 48h)
        hl_final = min(hl_val, 48.0)

        # Return the columns to be added
        return pd.Series([hl_final, source, proxy_cell], index=features)

    # 3. Apply to Main DataFrame
    # axis=1 passes the row to the function
    features_df = all_data.apply(_get_halflife_features, axis=1)

    # 4. Merge and Return
    # Concatenate the new columns to the main dataframe matching indices
    enriched_df = pd.concat([all_data, features_df], axis=1)

    print("Features populated successfully.")
    return enriched_df, features
