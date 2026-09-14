"""Figure 6: local fold stability and binding-site opening energy separate efficacy independently.

Both descriptors come from ViennaRNA. Local stability is the sliding-window MFE of the target
region; opening energy is the cost of holding a segment of the site single-stranded, computed
over nine segment lengths and two window anchorings. The two are correlated but not
interchangeable, so the figure splits ASOs into four groups at the median of each and compares
their within-experiment inhibition residual.

Prints the within-experiment Spearman of every fold and accessibility feature, the Pearson
matrix between the two families, and the four group means; writes the figure as a PNG.

    python notebooks/graphs/fold_vs_accessibility.py --matrix MATRIX.parquet
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

INHIBITION_PERCENT = "inhibition_percent"
CUSTOM_ID = "custom_id"

# The representatives the figure splits on: the strongest of each family by within-experiment
# Spearman. Every pair of one fold and one accessibility feature is correlated, from -0.48 to
# -0.82, so the figure names which pair it used.
FOLD_SPLIT = "fold_mfe_win40_flank30_step5"
ACCESS_SPLIT = "access_f60_sinf_u13_a5"

GROUPS = ["open / open", "open fold / buried site", "buried fold / open site", "buried / buried"]
COLORS = ["#2a7f62", "#8fbf9f", "#d8a373", "#b4573a"]


def feature_columns(matrix_path):
    """The fold and accessibility columns present in the matrix, in increasing scale order."""
    names = pq.read_schema(matrix_path).names
    fold = sorted(
        (c for c in names if c.startswith("fold_mfe_win")),
        key=lambda c: int(c.split("win")[1].split("_")[0]),
    )
    access = sorted(
        (c for c in names if c.startswith("access_") and c.endswith(("_a3", "_a5"))),
        key=lambda c: (int(c.split("_u")[1].split("_")[0]), c),
    )
    return fold, access


def within_experiment_spearman(df, column, min_rows=10):
    """Median Spearman of one feature with inhibition inside an experiment."""
    sub = df[[CUSTOM_ID, column, INHIBITION_PERCENT]].dropna()
    sizes = sub[CUSTOM_ID].map(sub[CUSTOM_ID].value_counts())
    sub = sub[sizes >= min_rows]
    if sub.empty:
        return float("nan"), 0

    experiment = sub[CUSTOM_ID]
    feature_rank = sub.groupby(experiment)[column].rank()
    inhibition_rank = sub.groupby(experiment)[INHIBITION_PERCENT].rank()
    feature_centred = feature_rank - feature_rank.groupby(experiment).transform("mean")
    inhibition_centred = inhibition_rank - inhibition_rank.groupby(experiment).transform("mean")

    covariance = (feature_centred * inhibition_centred).groupby(experiment).sum()
    spread = np.sqrt(
        (feature_centred**2).groupby(experiment).sum() * (inhibition_centred**2).groupby(experiment).sum()
    )
    # A group where either side is constant has no spread and no correlation to report.
    rho = (covariance / spread).replace([np.inf, -np.inf], np.nan).dropna()
    return float(rho.median()), len(rho)


def assign_groups(df, fold_column, access_column, within_experiment):
    """Label each row by whether the site is open or buried on each descriptor.

    A higher MFE is a less stable local fold and a lower opening energy is a more accessible
    site, so "open" is above the median on one and below it on the other.
    """
    if within_experiment:
        fold_cut = df.groupby(CUSTOM_ID)[fold_column].transform("median")
        access_cut = df.groupby(CUSTOM_ID)[access_column].transform("median")
    else:
        fold_cut = df[fold_column].median()
        access_cut = df[access_column].median()
    fold_open = df[fold_column] > fold_cut
    access_open = df[access_column] < access_cut
    return np.where(
        fold_open & access_open,
        GROUPS[0],
        np.where(
            fold_open & ~access_open,
            GROUPS[1],
            np.where(~fold_open & access_open, GROUPS[2], GROUPS[3]),
        ),
    )


def group_stats(df, labels):
    """Mean within-experiment residual and its 95% confidence interval, per group."""
    rows = []
    for name in GROUPS:
        residual = df.loc[labels == name, "residual"].to_numpy()
        mean = float(residual.mean())
        half_width = 1.96 * float(residual.std(ddof=1)) / np.sqrt(len(residual))
        rows.append({"group": name, "n": len(residual), "mean": mean, "ci": half_width})
    return rows


def draw(rows, fold_column, access_column, pearson, out_path):
    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    x = np.arange(len(rows))
    means = np.array([r["mean"] for r in rows])
    errors = np.array([r["ci"] for r in rows])

    low = float(min(means - errors))
    high = float(max(means + errors))
    span = high - low
    ax.set_ylim(low - 0.30 * span, high + 0.18 * span)
    bottom, top = ax.get_ylim()
    pad = 0.045 * (top - bottom)

    ax.axhline(0.0, color="#555555", linestyle="--", linewidth=1.0, zorder=1)
    ax.bar(x, means, yerr=errors, capsize=4, color=COLORS, edgecolor="#333333", linewidth=0.6, zorder=2)
    for xi, mean, error, row in zip(x, means, errors, rows):
        above = mean >= 0
        ax.text(
            xi,
            mean + error + pad if above else mean - error - pad,
            f"{mean:+.1f}",
            ha="center",
            va="bottom" if above else "top",
            fontsize=10,
            fontweight="bold",
        )
        ax.text(xi, bottom + pad, f"n={row['n']:,}", ha="center", va="bottom", fontsize=8, color="#555555")

    ax.set_xticks(x)
    ax.set_xticklabels([r["group"].replace(" / ", "\n") for r in rows], fontsize=9)
    ax.set_ylabel("within-experiment inhibition residual (pp)")
    ax.set_title("Fold stability and opening energy are complementary", fontsize=12)
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout(rect=(0, 0.075, 1, 1))
    fig.text(
        0.01,
        0.015,
        f"median split on {fold_column} and {access_column};  Pearson = {pearson:.2f};  "
        f"bars are group means \u00b1 95% CI",
        fontsize=7,
        color="#555555",
    )
    fig.savefig(out_path, dpi=300)
    print(f"\nwrote {out_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix", required=True, help="feature matrix parquet")
    ap.add_argument("--fold", default=FOLD_SPLIT, help="fold feature the figure splits on")
    ap.add_argument("--access", default=ACCESS_SPLIT, help="accessibility feature the figure splits on")
    ap.add_argument(
        "--within-experiment-split",
        action="store_true",
        help="take each median inside an experiment rather than over the whole dataset",
    )
    ap.add_argument(
        "--out",
        default=str(Path(__file__).with_name("fold_vs_accessibility.png")),
        help="where to write the figure",
    )
    args = ap.parse_args()

    fold, access = feature_columns(args.matrix)
    columns = list(dict.fromkeys(fold + access + [CUSTOM_ID, INHIBITION_PERCENT]))
    df = pd.read_parquet(args.matrix, columns=columns)
    df = df[df[INHIBITION_PERCENT].notna()].reset_index(drop=True)
    df["residual"] = df[INHIBITION_PERCENT] - df.groupby(CUSTOM_ID)[INHIBITION_PERCENT].transform("mean")
    print(f"rows={len(df):,}  experiments={df[CUSTOM_ID].nunique():,}")

    print("\nwithin-experiment Spearman with inhibition")
    for column in fold + access:
        rho, n = within_experiment_spearman(df, column)
        print(f"  {column:34s} {rho:+.4f}  ({n} experiments)")

    print("\nPearson between local stability and opening energy")
    correlations = df[fold + access].corr().loc[fold, access]
    header = "".join(f"{c.split('_u')[1]:>9s}" for c in access)
    print(f"  {'':34s}{header}")
    for fold_column in fold:
        cells = "".join(f"{correlations.at[fold_column, c]:9.3f}" for c in access)
        print(f"  {fold_column:34s}{cells}")

    split = df[[args.fold, args.access, CUSTOM_ID, "residual"]].dropna().reset_index(drop=True)
    pearson = float(np.corrcoef(split[args.fold], split[args.access])[0, 1])
    labels = assign_groups(split, args.fold, args.access, args.within_experiment_split)
    rows = group_stats(split, labels)

    print(f"\nfour-group split on {args.fold} x {args.access}")
    print(f"  Pearson = {pearson:.3f}")
    for row in rows:
        print(
            f"  {row['group']:26s} {row['n']:7d} {row['mean']:+7.2f}  "
            f"[{row['mean'] - row['ci']:+.2f}, {row['mean'] + row['ci']:+.2f}]"
        )

    draw(rows, args.fold, args.access, pearson, args.out)


if __name__ == "__main__":
    main()
