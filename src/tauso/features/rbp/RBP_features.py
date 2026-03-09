# ================= RNA-mode RBP_compatibility =================
# Works entirely in RNA alphabet (A/C/G/U). No U->T conversions.
import numpy as np
from scipy.stats import entropy


def create_positional_sequence_columns(
    df, source_col="flank_sequence_50", flank_size=50
):
    """
    Dynamically slices the sequence into Left, Core, and Right based on flank size.
    Left  = First [flank_size] nucleotides
    Right = Last [flank_size] nucleotides
    Core  = Everything in between (The ASO site)
    """
    seqs = df[source_col].fillna("").astype(str)

    # Slice
    df["seq_region_left"] = seqs.str[:flank_size]
    df["seq_region_right"] = seqs.str[-flank_size:]
    df["seq_region_core"] = seqs.str[flank_size:-flank_size]

    return df


def add_global_complexity_features(df, feature_cols, suffix):
    """Calculates Total Sum and Entropy for a list of columns."""
    matrix = df[feature_cols].fillna(0.0).values

    # 1. Total Sum
    total_col = f"RBP_interaction_total_{suffix}"
    total_scores = matrix.sum(axis=1)
    df[total_col] = total_scores

    # 2. Entropy
    with np.errstate(divide="ignore", invalid="ignore"):
        probs = matrix / total_scores[:, None]
        probs = np.nan_to_num(probs)  # Replace NaNs with 0

    div_col = f"RBP_diversity_global_{suffix}"
    df[div_col] = entropy(probs, axis=1)

    return df, [total_col, div_col]


def add_strict_functional_features(df, rbp_map, suffix):
    """Calculates Stabilizer/Destabilizer Sums and Contrasts."""
    new_feats = []
    role_cols = {}

    for role in ["stabilizer", "destabilizer"]:
        genes = [g for g, r in rbp_map.items() if r == role]
        # Reconstruct column names
        cols = [
            f"RBP_{g}_expr_aff_{suffix}"
            for g in genes
            if f"RBP_{g}_expr_aff_{suffix}" in df.columns
        ]

        if not cols:
            continue

        # Sum
        name = f"RBP_interaction_{role}_{suffix}"
        df[name] = df[cols].sum(axis=1)
        new_feats.append(name)
        role_cols[role] = name

    # Contrast
    if role_cols.get("stabilizer") and role_cols.get("destabilizer"):
        stab = df[role_cols["stabilizer"]]
        destab = df[role_cols["destabilizer"]]

        # Delta
        delta = f"RBP_interaction_delta_stab_vs_destab_{suffix}"
        df[delta] = stab - destab
        new_feats.append(delta)

        # Ratio
        ratio = f"RBP_interaction_ratio_stab_vs_destab_{suffix}"
        df[ratio] = stab / (destab + 1e-6)
        new_feats.append(ratio)

    return df, new_feats


def get_background_probs(target_sequence):
    """Calculates nt frequencies once."""
    if not target_sequence:
        return np.array([0.25, 0.25, 0.25, 0.25])

    # Fast count using string methods (faster than Counter for simple ACGT)
    s = target_sequence.upper()
    l = len(s) + 4  # +4 for pseudocounts

    p_a = (s.count("A") + 1) / l
    p_c = (s.count("C") + 1) / l
    p_g = (s.count("G") + 1) / l
    p_t = (s.count("T") + s.count("U") + 1) / l

    return np.array([p_a, p_c, p_g, p_t])
