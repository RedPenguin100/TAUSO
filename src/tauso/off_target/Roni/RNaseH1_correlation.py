import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import pearsonr, spearmanr
'''
# --- SETTINGS ---
# inhibition data
may_df = pd.read_csv("/home/oni/ASOdesign/scripts/data_genertion/data_updated_inhibition.csv")
inhibition = may_df[['index', 'Inhibition(%)']]

# RNaseH1 sizes to process
sizes = [0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800]

plt.figure(figsize=(10, 6))
results = []

for size in sizes:
    file = f"off_target_RNaseH1_{size}.csv"

    try:
        df = pd.read_csv(file)
    except Exception as e:
        print(f"❌ Could not read {file}: {e}")
        continue

    if "off_target_score" not in df.columns:
        print(f"⚠️ Missing off_target_score column in {file}")
        continue

    df = df[['index', 'off_target_score']]
    merged = pd.merge(inhibition, df, on="index", how="inner")

    # correlation including zeros
    if len(merged) > 1:
        corr_all, p_all = spearmanr(merged['Inhibition(%)'], merged['off_target_score'])
        results.append({"size": size, "mode": "all", "correlation": corr_all, "p_value": p_all})
        plt.scatter(size, corr_all,
                    color="blue",
                    marker='o',
                    label="all" if "all" not in plt.gca().get_legend_handles_labels()[1] else "")

    # correlation excluding zeros
    merged_nz = merged[merged['off_target_score'] != 0]
    if len(merged_nz) > 1:
        corr_nz, p_nz = spearmanr(merged_nz['Inhibition(%)'], merged_nz['off_target_score'])
        results.append({"size": size, "mode": "nonzero", "correlation": corr_nz, "p_value": p_nz})
        plt.scatter(size, corr_nz,
                    color="red",
                    marker='x',
                    label="nonzero" if "nonzero" not in plt.gca().get_legend_handles_labels()[1] else "")

# legend for marker meaning
handles, labels = plt.gca().get_legend_handles_labels()
custom_handles = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', label="all", markersize=8),
    Line2D([0], [0], marker='x', color='w', markeredgecolor='red', label="nonzero", markersize=8)
]

plt.xlabel("RNaseH1 size")
plt.ylabel("Correlation with Inhibition(%)")
plt.title("Spearman Correlation of RNaseH1 Off-target Scores with Inhibition")
plt.legend(handles=custom_handles, title="Mode", loc="best")
plt.grid(True)
plt.show()

# Create summary table
results_df = pd.DataFrame(results)
results_df = results_df.sort_values(by=["mode", "size"])
results_df["correlation"] = results_df["correlation"].round(5)

print(results_df)
results_df.to_csv("corr_pval_RNaseH1_spear.csv", index=False)
'''


# --- SETTINGS ---
inhibition = pd.read_csv("/home/oni/ASOdesign/scripts/data_genertion/data_updated_inhibition.csv")
df = pd.read_csv("../data_genertion/off_target/RNaseH/off_target_RNaseH1_0.csv")

# keep only relevant columns
df = df[['index', 'off_target_score']]
merged = pd.merge(inhibition, df, on="index", how="inner")

# group by Modification + true_length_of_seq
grouped = merged.groupby(["Modification", "true_length_of_seq"])

results = []
for (mod, length), group in grouped:
    if len(group) > 1:  # need at least 2 points for correlation
        corr, pval = spearmanr(group["Inhibition(%)"], group["off_target_score"])
        results.append({
            "group": f"{mod}/deoxy_{length}",
            "correlation": corr,
            "p_value": pval,
            "n_rows": len(group)   # keep group size
        })

results_df = pd.DataFrame(results)

# sort by correlation (optional)
results_df = results_df.sort_values(by="correlation", ascending=False)

# coloring: light red for positive, light blue for negative
colors = results_df["correlation"].apply(
    lambda x: "lightcoral" if x > 0 else "lightblue"
)

# --- PLOTTING ---
plt.figure(figsize=(12, 6))
plt.bar(results_df["group"], results_df["correlation"], color=colors, edgecolor="black")

plt.axhline(0, color="black", linewidth=0.8)
plt.xticks(rotation=90, fontsize=9)
plt.ylabel("Spearman Correlation with Inhibition(%)")
plt.title("Correlation by Modification × Length (RNaseH1 size 0)")

plt.tight_layout()
plt.grid(True)
plt.show()

# save table
results_df.to_csv("corr_by_mod_length_RNaseH1_0.csv", index=False)
print(results_df)