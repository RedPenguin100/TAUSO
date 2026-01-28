import pandas as pd

df = pd.read_csv("IDO1_exp_data.csv")

# -----------------------------
# 1. Check sequence consistency
# -----------------------------
seq_check = df.groupby("ASO")["Sequence"].nunique()
bad_asos = seq_check[seq_check > 1]

if not bad_asos.empty:
    raise ValueError(f"Inconsistent sequences found for ASO IDs: {bad_asos.index.tolist()}")

# -----------------------------
# 2. Add replicate counter
# -----------------------------
df["rep"] = df.groupby(["ASO", "Experiment"]).cumcount() + 1

# -----------------------------
# 3. Pivot to wide format
# -----------------------------
pivot_df = df.pivot(
    index=["ASO", "Sequence"],
    columns=["Experiment", "rep"],
    values="Inhibition rate"
)

# -----------------------------
# 4. Flatten column names
# -----------------------------
pivot_df.columns = [
    f"inhibition_rate_{exp.lower()}_{rep}" for exp, rep in pivot_df.columns
]

# -----------------------------
# 5. Reset index
# -----------------------------
pivot_df = pivot_df.reset_index()

# -----------------------------
# 6. Compute average inhibition rate
# -----------------------------
pivot_df["average_inhibition_rate"] = pivot_df.filter(like="inhibition_rate").mean(axis=1, skipna=True)


pivot_df = pivot_df.sort_values(by="average_inhibition_rate", ascending=False)

# ------------------------------
# 7. sequence length
# ------------------------------

pivot_df.insert(loc=pivot_df.columns.get_loc("Sequence")+1, column="Length", value=pivot_df["Sequence"].str.len())

# Save or inspect
print(pivot_df.head())
pivot_df.to_csv("IDO1_exp_inhibition.csv", index=False)