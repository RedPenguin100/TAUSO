"""
Process raw indexed oligo data into the filtered and averaged dataset used for model training.

Steps:
  1. Load indexed CSV
  2. Rename/standardize columns
  3. Compute structure features (sense_start, sense_length, etc.)
  4. Filter: remove unsupported chemistry, unmapped sequences, non-standard gapmers
  5. Save processed (pre-averaging) CSV
  6. Average duplicate measurements (same ASO, same experimental conditions)
  7. Save averaged CSV  ← this is what train_model.py reads
  8. Save diagnostic plots and tables to --output-dir

Usage:
    python -m notebooks.data.OligoAI.process_data --cpus 32
    python -m notebooks.data.OligoAI.process_data --cpus 8 --output-dir /path/to/plots

SLURM note: run with `python -u` or PYTHONUNBUFFERED=1 for live log output.
"""
import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
import pandas as pd

matplotlib.use("Agg")  # non-interactive backend for SLURM/headless
import matplotlib.pyplot as plt
import seaborn as sns

from notebooks.consts import (
    OLIGO_CSV_INDEXED,
    OLIGO_CSV_PROCESSED,
    OLIGO_CSV_PROCESSED_AVERAGED,
)
from notebooks.preprocessing import process_oligo_data, process_oligo_data_rename
from tauso.data.consts import INHIBITION, SEQUENCE
from tauso.populate.calculators.calculator import Calculator

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

def _plot_sequence_split_features(averaged_df, output_dir):
    """Bar chart: which features cause identical sequences to get multiple rows."""
    duplicates = averaged_df[averaged_df.duplicated(subset=[SEQUENCE], keep=False)]
    variation_counts = defaultdict(int)

    for _, group in duplicates.groupby(SEQUENCE):
        for col in group.columns:
            if col in [SEQUENCE, INHIBITION, "index_oligo"]:
                continue
            if group[col].nunique(dropna=False) > 1:
                variation_counts[col] += 1

    variation_df = pd.DataFrame(
        list(variation_counts.items()), columns=["Feature", "Number of Sequences Split"]
    ).sort_values(by="Number of Sequences Split", ascending=False)

    variation_df.to_csv(output_dir / "sequence_split_features.csv", index=False)
    logger.info("Saved sequence split table -> %s", output_dir / "sequence_split_features.csv")

    if variation_df.empty:
        logger.info("No identical sequences split by any feature — skipping plot.")
        return

    plt.figure(figsize=(10, max(4, len(variation_df) * 0.4)))
    sns.barplot(data=variation_df, x="Number of Sequences Split", y="Feature", palette="magma")
    plt.title("Which Features Cause Identical Sequences to Split?")
    plt.xlabel("Number of Unique Sequences Affected")
    plt.ylabel("Varying Feature")
    plt.tight_layout()
    plt.savefig(output_dir / "sequence_variation_plot.png", dpi=150)
    plt.close()
    logger.info("Saved sequence variation plot -> %s", output_dir / "sequence_variation_plot.png")


def _analyze_sequence_custom_id_splits(averaged_df, output_dir):
    """Report which columns cause rows with identical (SEQUENCE, custom_id) to differ."""
    subset_cols = [SEQUENCE, "custom_id"]
    duplicates = averaged_df[averaged_df.duplicated(subset=subset_cols, keep=False)]

    if duplicates.empty:
        logger.info("No (SEQUENCE, custom_id) duplicates found.")
        return

    grouped = duplicates.groupby(subset_cols)
    total = len(grouped)
    variation_counts = defaultdict(int)

    for _, group in grouped:
        for col in group.columns:
            if col in subset_cols + [INHIBITION, "index_oligo"]:
                continue
            if group[col].nunique(dropna=False) > 1:
                variation_counts[col] += 1

    rows = [
        {"Feature": col, "Pairs_Split": cnt, "Pct": cnt / total * 100}
        for col, cnt in sorted(variation_counts.items(), key=lambda x: x[1], reverse=True)
    ]
    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(output_dir / "custom_id_split_features.csv", index=False)

    logger.info("(SEQUENCE, custom_id) duplicate pairs: %d", total)
    for r in rows:
        logger.info("  %-25s | split %d pairs (%.1f%%)", r["Feature"], r["Pairs_Split"], r["Pct"])
    logger.info(
        "Saved custom_id split table -> %s", output_dir / "custom_id_split_features.csv"
    )


def _save_cell_line_summary(processed_df, output_dir):
    counts = processed_df["Cell_line"].value_counts().reset_index()
    counts.columns = ["Cell_line", "Count"]
    counts.to_csv(output_dir / "cell_line_counts.csv", index=False)
    logger.info("Cell lines present (%d):", len(counts))
    for _, row in counts.iterrows():
        logger.info("  %-35s %d", row["Cell_line"], row["Count"])
    logger.info("Saved cell line summary -> %s", output_dir / "cell_line_counts.csv")


def _plot_inhibition_distribution(averaged_df, output_dir):
    plt.figure(figsize=(8, 4))
    plt.hist(averaged_df[INHIBITION].dropna(), bins=80, edgecolor="none", color="steelblue")
    plt.xlabel("Inhibition (%)")
    plt.ylabel("Count")
    plt.title(f"Inhibition distribution (n={len(averaged_df):,})")
    plt.tight_layout()
    plt.savefig(output_dir / "inhibition_distribution.png", dpi=150)
    plt.close()
    logger.info("Saved inhibition distribution -> %s", output_dir / "inhibition_distribution.png")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--cpus", type=int, default=1, help="CPUs for structure calculation")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).parent / "process_data_output",
        help="Directory for plots and diagnostic tables",
    )
    parser.add_argument(
        "--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"]
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s | %(message)s",
        stream=sys.stdout,
        force=True,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Load
    # ------------------------------------------------------------------
    logger.info("Loading indexed data from %s", OLIGO_CSV_INDEXED)
    data = pd.read_csv(OLIGO_CSV_INDEXED)
    logger.info("Loaded %d rows", len(data))

    # ------------------------------------------------------------------
    # 2. Rename / standardize columns
    # ------------------------------------------------------------------
    logger.info("Renaming columns...")
    data = process_oligo_data_rename(data)

    # ------------------------------------------------------------------
    # 3. Structure features (sense_start, sense_length, etc.)
    # ------------------------------------------------------------------
    logger.info("Computing structure features with %d CPUs...", args.cpus)
    calculator = Calculator(data=data, data_version="oligo", overwrite=True, cpus=args.cpus)
    calculator.calculate_structure()
    data = calculator.data

    # ------------------------------------------------------------------
    # 4. Filter
    # ------------------------------------------------------------------
    logger.info("Filtering (strict gapmer patterns)...")
    processed_data = process_oligo_data(data, strict_gapmer_patterns=True)
    logger.info("After filtering: %d rows", len(processed_data))

    # ------------------------------------------------------------------
    # 5. Save processed (pre-averaging)
    # ------------------------------------------------------------------
    logger.info("Saving processed data -> %s", OLIGO_CSV_PROCESSED)
    processed_data.to_csv(OLIGO_CSV_PROCESSED)

    _save_cell_line_summary(processed_data, args.output_dir)

    # ------------------------------------------------------------------
    # 6. Average duplicate measurements
    # ------------------------------------------------------------------
    logger.info("Averaging duplicate measurements...")
    feature_cols = [col for col in processed_data.columns if col not in [INHIBITION, "index_oligo"]]
    averaged_df = processed_data.groupby(feature_cols, as_index=False, dropna=False).agg(
        {INHIBITION: "mean", "index_oligo": "first"}
    )
    final_cols = ["index_oligo"] + feature_cols + [INHIBITION]
    averaged_df = averaged_df[final_cols]
    logger.info("Collapsed %d -> %d rows after averaging", len(processed_data), len(averaged_df))

    # ------------------------------------------------------------------
    # 7. Save averaged
    # ------------------------------------------------------------------
    logger.info("Saving averaged data -> %s", OLIGO_CSV_PROCESSED_AVERAGED)
    averaged_df.to_csv(OLIGO_CSV_PROCESSED_AVERAGED, index=False)

    # ------------------------------------------------------------------
    # 8. Diagnostics
    # ------------------------------------------------------------------
    logger.info("Running diagnostics...")
    _plot_sequence_split_features(averaged_df, args.output_dir)
    _analyze_sequence_custom_id_splits(averaged_df, args.output_dir)
    _plot_inhibition_distribution(averaged_df, args.output_dir)

    logger.info("Done. Outputs written to %s", args.output_dir)


if __name__ == "__main__":
    main()
