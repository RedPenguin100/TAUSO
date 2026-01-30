import os
import re
from pathlib import Path

import pandas as pd
from pandarallel import pandarallel

from ...common.gtf import filter_gtf_genes
from ...data.data import load_db, get_paths
from ...file_utils import get_fasta_dict_from_path
from ...genome.read_human_genome import get_locus_to_data_dict


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


def load_cell_line_gene_maps(cell_map, data_dir, valid_db_genes,
                             n_specific=300, n_fallback_scan=500,
                             filter_mode='protein_coding',
                             genome_db=None):
    """
    Loads expression data, filtering strictly against the GTF/GFF database annotations.
    """

    # 1. Initialize DB Connection
    if genome_db is None:
        # Calls your existing function
        genome_db = load_db()

        # 2. Build the Biological Allowlist (The "Truth" Set)
    gtf_allowed_genes = filter_gtf_genes(genome_db, filter_mode)

    # 3. Intersect with CoordinateMapper genes
    # The gene must be biologically valid (GTF) AND technically valid (in your DB/Mapper)
    final_valid_genes = set(valid_db_genes).intersection(gtf_allowed_genes)

    print(f"Filter Ready: {len(final_valid_genes)} genes passed specificiation ({filter_mode}).")
    print("Processing cell lines...")

    cell_line_top_genes = {}
    all_scan_sets = []
    global_fetch_set = set()

    # Regex to clean symbols like "GAPDH (2597)" -> "GAPDH"
    name_fix_regex = re.compile(r" \(\d+\)")
    n_scan = max(n_specific, n_fallback_scan)

    for cell_name, ach_id in cell_map.items():
        p = data_dir / f"{ach_id}_expression.csv"
        if not p.exists():
            continue

        # Load and sort
        df = pd.read_csv(p).sort_values('expression_TPM', ascending=False)

        # A. Clean Names
        df['Gene'] = df['Gene'].str.replace(name_fix_regex, '', regex=True).str.strip()

        # B. Apply The "Truth" Filter
        # This one step replaces all prefix checks, pseudogene regexes, and exclusion lists.
        df = df[df['Gene'].isin(final_valid_genes)]

        # C. Fetch Scan Set
        genes_scan = df.head(n_scan)['Gene'].tolist()

        if not genes_scan:
            print(f"  -> Warning: No valid genes found for {cell_name}")
            continue

        # 1. Specific Set
        genes_specific = genes_scan[:n_specific]
        cell_line_top_genes[cell_name] = genes_specific
        global_fetch_set.update(genes_specific)

        # 2. Intersection Set
        all_scan_sets.append(set(genes_scan[:n_fallback_scan]))

        print(f"  -> {cell_name}: Loaded top {len(genes_specific)} genes.")

    # --- Generate Consensus ---
    fallback_genes = []
    if all_scan_sets:
        consensus_genes = set.intersection(*all_scan_sets)
        fallback_genes = list(consensus_genes)
        print(f"\nFallback Consensus Intersection size: {len(fallback_genes)}")
        global_fetch_set.update(fallback_genes)
    else:
        print("\nWarning: No data available for intersection.")

    return cell_line_top_genes, fallback_genes, global_fetch_set


def load_cell_line_transcriptomes(cell_map, db, expression_dir, filter_mode='non_mt',
                                  add_transcript=False):
    valid_genes = filter_gtf_genes(db, filter_mode=filter_mode)

    # Lazy Load: Only load the heavy sequence map if we actually need it
    locus_map = None
    if add_transcript:
        locus_map = get_locus_to_data_dict(include_introns=True)

    transcriptomes = {}

    for cell_name, ach_id in cell_map.items():
        path = os.path.join(expression_dir, f"{ach_id}_expression.csv")
        if not os.path.exists(path):
            print(f"  -> Warning: No expression data for {cell_name}")
            continue

        # Load and Filter
        df = pd.read_csv(path)
        df = df[df['Gene'].isin(valid_genes)].copy()

        # Transcript Logic (Conditional)
        if add_transcript:
            df["Original Transcript Sequence"] = df["Gene"].map(
                lambda g: locus_map[g].full_mrna.replace('T', 'U') if g in locus_map else None
            )
            df = df.dropna(subset=["Original Transcript Sequence"])

        if not df.empty:
            transcriptomes[cell_name] = df

    return transcriptomes