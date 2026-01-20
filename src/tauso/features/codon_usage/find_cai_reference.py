import re

import pandas as pd
from pathlib import Path


# --- CORE FUNCTIONS ---

def get_top_expressed_genes(file_path, n_top=300):
    """
    Reads an expression CSV, filters out mitochondrial genes,
    and returns the top N most highly expressed nuclear genes.
    """
    if not file_path.exists():
        return []

    try:
        # 1. Smart load: scan header for the 'expression_norm' column
        header = pd.read_csv(file_path, nrows=0).columns.tolist()
        target_col = next((c for c in header if "expression_norm" in c), None)

        if not target_col:
            print(f"Warning: 'expression_norm' column missing in {file_path.name}")
            return []

        # 2. Load only necessary columns
        df = pd.read_csv(file_path, usecols=['Gene', target_col])

        # 3. Filter Nuclear (Exclude 'MT-') & Sort
        # Ensure we drop NaNs in the target column before sorting
        df = df.dropna(subset=[target_col])
        nuclear_df = df[~df['Gene'].str.startswith('MT-')]

        top_genes = (
            nuclear_df.sort_values(target_col, ascending=False)
            .head(n_top)['Gene']
            .apply(lambda x: x.split(' ')[0])  # Clean "ACTB (60)" -> "ACTB"
            .tolist()
        )
        return top_genes

    except Exception as e:
        print(f"Error reading {file_path.name}: {e}")
        return []


def load_cell_line_gene_maps(cell_map, data_dir, valid_db_genes, n_specific=300, n_fallback_scan=500):
    """
    Fixed version: Filters for nuclear protein-coding genes present in the DB
    before calculating top genes and intersections.
    """
    print("Identifying top valid protein-coding genes per cell line...")

    cell_line_top_genes = {}
    all_scan_sets = []
    global_fetch_set = set()

    # Regex to clean symbols like "GAPDH (2597)" -> "GAPDH"
    name_fix_regex = re.compile(r" \(\d+\)")

    # Biological exclusion filters
    non_coding_prefixes = ('MT-', 'MIR', 'SNO', 'RNU', 'LINC', 'SNOR', 'ENSG', 'RN7', 'RNA5')
    pseudogene_pattern = re.compile(r'-PS\d+$|\d+P\d+$|P\d+$')

    # Ensure scan size covers both needs
    n_scan = max(n_specific, n_fallback_scan)

    for cell_name, ach_id in cell_map.items():
        p = data_dir / f"{ach_id}_expression.csv"
        if not p.exists():
            continue

        # Load and sort once
        df = pd.read_csv(p).sort_values('expression_TPM', ascending=False)

        # A. Clean Names
        df['Gene'] = df['Gene'].str.replace(name_fix_regex, '', regex=True).str.strip()

        # B. Biological Filters
        df = df[~df['Gene'].str.startswith(non_coding_prefixes, na=False)]
        df = df[~df['Gene'].str.contains(pseudogene_pattern, regex=True, na=False)]

        # C. Database Cross-Reference
        # Only keep genes that the CoordinateMapper knows about
        df = df[df['Gene'].isin(valid_db_genes)]

        # D. Fetch the Top Scan Set (vetted genes)
        genes_scan = df.head(n_scan)['Gene'].tolist()

        if not genes_scan:
            print(f"  -> Warning: No valid genes found for {cell_name}")
            continue

        # 1. Specific Set (Top 300)
        genes_specific = genes_scan[:n_specific]
        cell_line_top_genes[cell_name] = genes_specific
        global_fetch_set.update(genes_specific)

        # 2. Store for Intersection (Top 500)
        all_scan_sets.append(set(genes_scan[:n_fallback_scan]))

        print(f"  -> {cell_name}: Loaded top {len(genes_specific)} specific genes.")

    # --- Generate Consensus Fallback (Intersection) ---
    fallback_genes = []
    if all_scan_sets:
        consensus_genes = set.intersection(*all_scan_sets)

        # We take a list from the intersection
        fallback_genes = list(consensus_genes)
        print(f"\nFallback Consensus Intersection size: {len(fallback_genes)}")

        # Add fallback genes to global fetch queue
        global_fetch_set.update(fallback_genes)
    else:
        print("\nWarning: No data available for intersection.")

    return cell_line_top_genes, fallback_genes, global_fetch_set