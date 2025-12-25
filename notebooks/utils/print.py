from scipy.stats import pearsonr, spearmanr
from sklearn.feature_selection import mutual_info_regression


def print_correlations(df, name1, name2, p_value_threshold=None):
    if p_value_threshold is None:
        p_value_threshold = 1
    corr, p_value = pearsonr(df[name1], df[name2])
    if p_value < p_value_threshold:
        print(f"Feature: {name1:<35}, Pearson: {corr:<5.2f}, p-value: {p_value:<10.2} Target: {name2:<35}")
    corr, p_value = spearmanr(df[name1], df[name2])
    if p_value < p_value_threshold:
        print(f"Feature: {name1:<35}, Spearman: {corr:<5.2f}, p-value: {p_value:<10.2} Target: {name2:<35}")

def calc_mutual_information(df, x_col, y_col):
    """
    Calculates Mutual Information (MI) between two continuous variables.
    MI >= 0. Higher values indicate stronger dependency (linear or non-linear).
    """
    # 1. Drop rows where either column is NaN (MI cannot handle NaNs)
    clean_df = df[[x_col, y_col]].dropna()

    if len(clean_df) < 2:
        return 0.0

    # 2. Reshape X for sklearn (requires 2D array)
    X = clean_df[[x_col]]
    y = clean_df[y_col]

    # 3. Calculate MI
    # discrete_features=False tells it these are continuous numbers
    mi_score = mutual_info_regression(X, y, discrete_features=False, random_state=42)

    return mi_score[0]