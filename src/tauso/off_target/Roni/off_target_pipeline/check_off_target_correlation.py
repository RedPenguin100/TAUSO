import pandas as pd
import matplotlib.pyplot as plt
import glob
import re
from matplotlib.lines import Line2D
from scipy.stats import pearsonr

# --- SETTINGS ---
may_df = pd.read_csv("/home/oni/ASOdesign/scripts/data_genertion/data_updated_inhibition.csv")
inhibition = may_df[['index', 'Inhibition(%)']]

method = "TPM"

# file naming scheme: off_target.GEO.top{N}.cutoff{C}.premRNA.csv
#pattern = f"/home/oni/ASOdesign/scripts/data_genertion/off_target/off_target.{method}.top*.cutoff*.premRNA.csv"
pattern = f"/home/oni/ASOdesign/scripts/data_genertion/off_target/off_target.top*.cutoff*.premRNA.csv"

# color mapping
color_map = {
    "100": "blue",
    "200": "green",
    "500": "red",
}

plt.figure(figsize=(10, 6))

# store results for table
results = []

for file in glob.glob(pattern):
    # extract topN and cutoff from filename
    m = re.search(r'top(\d+)\.cutoff(\d+)', file)
    if not m:
        continue
    topN, cutoff = m.groups()

    df = pd.read_csv(file)
    df = df[['index', f'off_target_score_{method}']]
    merged = pd.merge(inhibition, df, on="index", how="inner")

    # correlation including zeros
    if len(merged) > 1:
        corr_all, p_all = pearsonr(merged['Inhibition(%)'], merged[f'off_target_score_{method}'])
        results.append({"topN": int(topN), "cutoff": int(cutoff),
                        "mode": f"{method}", "correlation": corr_all, "p_value": p_all})
        plt.scatter(int(cutoff), corr_all,
                    color=color_map.get(topN, "gray"),
                    marker='o',
                    label=f"top {topN}" if f"top {topN}" not in plt.gca().get_legend_handles_labels()[1] else "")

    # correlation excluding zeros
    merged_nz = merged[merged[f'off_target_score_{method}'] != 0]
    if len(merged_nz) > 1:
        corr_nz, p_nz = pearsonr(merged_nz['Inhibition(%)'], merged_nz[f'off_target_score_{method}'])
        results.append({"topN": int(topN), "cutoff": int(cutoff),
                        "mode": f"{method}_nonzero", "correlation": corr_nz, "p_value": p_nz})
        plt.scatter(int(cutoff), corr_nz,
                    color=color_map.get(topN, "gray"),
                    marker='x',
                    label="")

# add custom legend entries for marker meaning
handles, labels = plt.gca().get_legend_handles_labels()
custom_handles = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='black', label=f"{method}", markersize=8),
    Line2D([0], [0], marker='x', color='w', markeredgecolor='black', label=f"{method}_nonzero", markersize=8)
]

plt.xlabel("Cutoff")
plt.ylabel("Correlation with Inhibition(%)")
plt.title("Correlation of Off-target Scores with Inhibition (method - Arithmetic)")
plt.legend(handles=handles + custom_handles, title="Legend", loc="best")
plt.grid(True)
plt.show()

# Create summary table
# =========================
results_df = pd.DataFrame(results)

# enforce ordering: first GEO, then GEO_nonzero
results_df["mode"] = pd.Categorical(
    results_df["mode"],
    categories=[f"{method}", f"{method}_nonzero"],
    ordered=True
)

results_df = results_df.sort_values(by=["mode", "topN", "cutoff"])
results_df["correlation"] = results_df["correlation"].round(5)
print(results_df)
results_df.to_csv("corr_pval.csv", index=False)
