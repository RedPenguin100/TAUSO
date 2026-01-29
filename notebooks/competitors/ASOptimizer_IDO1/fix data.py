import pandas as pd

# Load the two datasets
aso_155 = pd.read_csv("/home/roniz/PycharmProjects/ASOdesign/scripts/Roni/ASOptimizer_IDO1/IDO1_ASO_155.csv")
model_scores = pd.read_csv("/home/roniz/PycharmProjects/ASOdesign/scripts/Roni/ASOptimizer_IDO1/IDO1_ASO_for_model_with_score.csv")

# In case the model_scores file has the first row as header accidentally:
if list(model_scores.columns)[0].startswith("C"):
    model_scores = pd.read_csv("/home/roniz/PycharmProjects/ASOdesign/scripts/Roni/ASOptimizer_IDO1/IDO1_ASO_for_model_with_score.csv", header=None)
    model_scores.columns = model_scores.iloc[0]
    model_scores = model_scores.drop(0).reset_index(drop=True)

# Convert column types for matching (ensure ASO IDs are strings)
aso_155["ASO"] = aso_155["ASO"].astype(str)
model_scores["ASO"] = model_scores["ASO"].astype(str)

# Merge score column by ASO
merged = pd.merge(aso_155, model_scores[["ASO", "score"]], on="ASO", how="left")

# Save result
merged.to_csv("IDO1_ASO_155_with_scores.csv", index=False)

print(f"✅ Merged successfully — total rows: {len(merged)}")
print(f"Matched scores: {merged['score'].notna().sum()} / {len(merged)}")
print("Saved to IDO1_ASO_155_with_scores.csv")
