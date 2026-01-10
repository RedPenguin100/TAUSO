import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from pathlib import Path

# Style parameters
sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12


def load_and_merge_feature(main_df, feature_name, features_dir):
    """
    Loads a specific feature CSV and merges it with the main dataset.

    Args:
        main_df (pd.DataFrame): The preprocessed main dataset (must have 'log_inhibition').
                                Assumes the index is the unique identifier (rna_id).
        feature_name (str): The name of the feature (and the CSV file name).
        features_dir (str or Path): Path to the folder containing feature CSVs.

    Returns:
        pd.DataFrame: Merged DataFrame containing only valid rows (intersection).
    """

    file_path = Path(features_dir) / f"{feature_name}.csv"
    if not file_path.exists():
        raise FileNotFoundError(f"Feature file not found: {file_path}")

    # Load feature data based on save_feature
    feature_df = pd.read_csv(file_path)

    # Ensure the 'index' column is set as the index for proper merging
    if 'index' not in feature_df.columns:
        raise ValueError(f"The file {feature_name}.csv is missing the 'index' column.")

    if 'index' not in main_df.columns:
        raise ValueError(f"Main df is missing the 'index' column.")

    # Merge with main_df (Inner join to keep only matching IDs)
    merged_df = main_df.merge(feature_df, on='index', how='inner')

    # Remove NaN values for the specific feature or inhibition
    merged_df = merged_df.dropna(subset=['log_inhibition', feature_name])

    return merged_df


def _add_stats_text(ax, x, y, x_label, y_label):
    """Helper function to calculate and display Spearman & Pearson stats."""

    # Pearson (Linear)
    r_p, p_p = stats.pearsonr(x, y)
    # Spearman (Rank / Monotonic)
    r_s, p_s = stats.spearmanr(x, y)

    # Formatting P-values
    def fmt_p(p):
        return "< 0.001" if p < 0.001 else f"{p:.3f}"

    text_str = (
        f"Pearson $r = {r_p:.2f}$ ($p$ {fmt_p(p_p)})\n"
        f"Spearman $\\rho = {r_s:.2f}$ ($p$ {fmt_p(p_s)})"
    )

    x_pos = 0.95 if r_p < 0 else 0.05
    ha = 'right' if r_p < 0 else 'left'
    props = dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray')
    ax.text(x_pos, 0.95, text_str, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', horizontalalignment=ha, bbox=props)


def plot_feature_scatter(df, feature_name,
                         x_label=None, y_label="Log Inhibition",
                         title=None, ax=None):
    """
    Creates a scatter plot with a regression line and statistical annotation.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))

    x = df[feature_name]
    y = df['log_inhibition']

    # Scatter plot with regression line
    sns.regplot(
        x=x, y=y, ax=ax,
        scatter_kws={'alpha': 0.5, 's': 15, 'edgecolor': 'none'},
        line_kws={'color': 'black', 'linewidth': 2}
    )

    # Add statistics
    _add_stats_text(ax, x, y, feature_name, 'log_inhibition')

    # Formatting
    ax.set_xlabel(x_label if x_label else feature_name)
    ax.set_ylabel(y_label)
    if title:
        ax.set_title(title)

    ax.grid(True, linestyle='--', alpha=0.7)

    return ax


def plot_binned_trend(df, feature_name, bins=10,x_label=None,
                      y_label="Mean Log Inhibition",title=None, ax=None):
    """
    Binds the feature values into quantiles and plots the mean inhibition per bin.
    Useful for visualizing trends in noisy data.
    """

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    local_df = df.copy()

    x = df[feature_name]
    y = df['log_inhibition']

    # 1. Binning (Using qcut for roughly equal number of points per bin)
    # If qcut fails (too many duplicate values), fall back to cut (equal width)

    try:
        local_df['bin'] = pd.qcut(local_df[feature_name], q=bins, duplicates='drop')
    except ValueError:
        local_df['bin'] = pd.cut(local_df[feature_name], bins=bins)

    # 2. Calculate means for each bin
    bin_stats = local_df.groupby('bin', observed=True).agg({
    feature_name: 'mean', 'log_inhibition': ['mean', 'std', 'count']}).reset_index()
    bin_stats.columns = ['bin', 'feature_mean', 'inhibition_mean', 'inhibition_std', 'count']
    bin_stats['inhibition_sem'] = bin_stats['inhibition_std'] / np.sqrt(bin_stats['count'])

    # Clean up (remove empty bins)
    bin_stats = bin_stats.dropna()
    X_binned = bin_stats['feature_mean']
    Y_binned = bin_stats['inhibition_mean']
    Y_error = bin_stats['inhibition_sem']

    # 3. Plot
    # Points representing the bins
    ax.errorbar(
        X_binned, Y_binned, yerr=Y_error,
        fmt='o', color='blue', ecolor='black', elinewidth=2, capsize=4,
        markersize=8, alpha=0.9, zorder=5,
        label='Bin Mean $\pm$ S.E.M.'
    )

    # Trend line through the bins
    sns.regplot(x=X_binned, y=Y_binned, ax=ax, scatter=False,
    line_kws={'color': 'red', 'linestyle': '--', 'linewidth': 2},
    label='Linear Fit (95% CI)'
    )

    # 4. Statistics (Calculated on the BINNED means, as requested)
    _add_stats_text(ax, x, y, feature_name, 'log_inhibition')

    # Formatting
    ax.set_xlabel(x_label if x_label else f"{feature_name} (Binned)")
    ax.set_ylabel(y_label)

    if title:
        ax.set_title(title)

    else:
        ax.set_title(f"Binned Trend: {feature_name} ({len(bin_stats)} bins)")

    ax.grid(True, linestyle='--', alpha=0.7)
    ax.legend(loc='best', fontsize=10)

    return ax

