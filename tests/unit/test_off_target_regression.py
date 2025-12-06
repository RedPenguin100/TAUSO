#!/usr/bin/env python3
"""
Script to run a real off-target search against the installed genome.
Usage: python scripts/run_real_search.py [SEQUENCE]
"""

import sys
import os
import argparse
import time
import pandas as pd

from tauso.off_target.search import find_all_gene_off_targets
from tauso.data import get_paths

def test_regression():
    # parser = argparse.ArgumentParser(description="Search for ASO off-targets in the real genome.")
    # # Default is the KLKB1 fragment
    # parser.add_argument("sequence", nargs="?", default="AGTGCCACATTAGAACAGCT",
    #                     help="ASO sequence to search (default: KLKB1 fragment)")
    # parser.add_argument("--mismatches", "-m", type=int, default=3, help="Max mismatches allowed")
    #
    # args = parser.parse_args()

    sequence = 'AGTGCCACATTAGAACAGCT'
    mismatches = 3

    # 1. Verify Environment
    print("Checking genome installation...")
    try:
        paths = get_paths()
        if not os.path.exists(paths['fasta']):
            print(f"Error: Genome FASTA not found at {paths['fasta']}")
            print("Please run 'tauso setup-genome' first.")
            sys.exit(1)
        print(f"Target Genome: {paths['fasta']}")
    except Exception as e:
        print(f"Error locating genome data: {e}")
        sys.exit(1)

    # 2. Run Search
    print(f"\n--- Starting Off-Target Search ---")
    print(f"Query: {sequence}")
    print(f"Length: {len(sequence)} bp")
    print(f"Max Mismatches: {mismatches}")
    print("Loading index and searching...")

    start_time = time.time()

    try:
        # Use the high-level function that returns a DataFrame
        hits_df = find_all_gene_off_targets(sequence, genome='GRCh38', max_mismatches=mismatches)
    except Exception as e:
        print(f"\nCRITICAL FAILURE during search: {e}")
        sys.exit(1)

    end_time = time.time()
    duration = end_time - start_time

    # 3. Report Results
    print(f"\n--- Results ---")
    print(f"Time taken: {duration:.4f} seconds")

    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 1000)
    print(hits_df)


    # FIX: Check if DataFrame is empty using .empty attribute
    if not hits_df.empty:
        print(f"Total Hits found: {len(hits_df)}")

        # Sort by mismatches (best first)
        hits_df = hits_df.sort_values(by=['mismatches', 'chrom', 'start'])

        print(f"\n{'Chrom':<10} | {'Start':<12} | {'End':<12} | {'Str':<3} | {'MM':<2} | {'Gene':<15} | {'Region'}")
        print("-" * 80)

        # Iterate over rows
        count = 0
        for _, hit in hits_df.iterrows():
            if count >= 20:
                print(f"... and {len(hits_df) - 20} more.")
                break

            gene = str(hit['gene_name']) if hit['gene_name'] else "None"
            print(
                f"{hit['chrom']:<10} | {hit['start']:<12} | {hit['end']:<12} | {hit['strand']:<3} | {hit['mismatches']:<2} | {gene:<15} | {hit['region_type']}")
            count += 1
    else:
        print("No matches found within mismatch threshold.")
