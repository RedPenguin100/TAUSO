"""Plot the per-cohort sequence-diversity ratio that drives the modification-scan filter.

Reuses the cleaning pipeline up to the chemistry filter, computes the per-cohort
diversity stats, and saves a horizontal bar chart annotated with the threshold and
the unique-design count per cohort.

Usage:
    python -m notebooks.data.ASOptimizer.plot_diversity
"""
import logging

import matplotlib

matplotlib.use("Agg")  # non-interactive backend for headless runs
import matplotlib.pyplot as plt
import seaborn as sns
from notebooks.data.ASOptimizer.consts import ASOPTIMIZER_DATA_PATH
from notebooks.data.ASOptimizer.data_cleanup import (
    MOD_SCAN_THRESHOLD,
    add_cohort,
    compute_cohort_diversity,
    run_pipeline,
)

logger = logging.getLogger(__name__)

OUTPUT_PNG = ASOPTIMIZER_DATA_PATH.parent / "diversity_ratio.png"


def plot_diversity(diversity_stats, threshold=MOD_SCAN_THRESHOLD, output_png=OUTPUT_PNG):
    plot_data = diversity_stats.sort_values("Diversity_Ratio")

    plt.figure(figsize=(12, 10))
    ax = sns.barplot(data=plot_data, x="Diversity_Ratio", y="Cohort", palette="viridis")
    plt.axvline(x=threshold, color="red", linestyle="--", label=f"Threshold ({threshold})")

    for i, patch in enumerate(ax.patches):
        count = plot_data.iloc[i]["Unique_Molecular_Designs"]
        ax.annotate(
            f" {int(count)} designs",
            (patch.get_width(), patch.get_y() + patch.get_height() / 2),
            va="center",
            fontsize=9,
            color="black",
            fontweight="bold",
        )

    plt.title("Sequence Diversity Ratio per Cohort (After Cross-Probe & Density Aggregation)")
    plt.xlabel("Diversity Ratio (Unique Sequences / Unique Molecular Designs)")
    plt.ylabel("Cohort (Gene + Cell Line)")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(output_png, dpi=150)
    plt.close()
    logger.info("saved %s", output_png)


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s | %(message)s")
    df = run_pipeline(through_mod_scan=False)
    df = add_cohort(df)
    diversity = compute_cohort_diversity(df)
    plot_diversity(diversity)


if __name__ == "__main__":
    main()
