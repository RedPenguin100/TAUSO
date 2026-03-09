import numpy as np

from ...algorithms.suffix_array import build_suffix_array, longest_prefix


def calculate_chimera_ars(suffix_array, target_sequence, step_size):
    longest_prefix_lengths = []

    for start_index in range(1, len(target_sequence), step_size):
        prefix, _ = longest_prefix(target_sequence[start_index:], suffix_array)
        longest_prefix_lengths.append(len(prefix))

    chimera_ars_score = np.mean(longest_prefix_lengths)
    return chimera_ars_score


def clean_region_from_aso(region, aso):
    if not isinstance(region, str) or not isinstance(aso, str):
        return region
    count = region.count(aso)
    return region.replace(aso, "") if count == 1 else region


def add_local_chimera_features_cds(
        df,
        flank_windows,
        sequence_col,
        local_col_template="local_coding_region_around_ASO_{flank}",
        out_col_template="chimera_score_CDS_{flank}",
        step_size=1,
        normalize_by_aso_len=False,
        suffix_array_cache=None,
):
    """
    Adds local CDS chimera scores for each flank window:
      local_coding_region_around_ASO_{flank} -> chimera_score_CDS_{flank}

    Returns (df, suffix_array_cache).
    """
    if suffix_array_cache is None:
        suffix_array_cache = {}

    def safe_build_suffix_array(region_clean: str):
        if region_clean in suffix_array_cache:
            return suffix_array_cache[region_clean]
        sa = build_suffix_array(region_clean)
        suffix_array_cache[region_clean] = sa
        return sa

    def compute_chimera(region, aso):
        if not isinstance(region, str) or not region.strip():
            return np.nan
        if not isinstance(aso, str) or not aso.strip():
            return np.nan

        region_clean = clean_region_from_aso(region, aso)
        if not isinstance(region_clean, str) or not region_clean.strip():
            return np.nan

        sa = safe_build_suffix_array(region_clean)
        score = calculate_chimera_ars(sa, aso, step_size=step_size)

        if normalize_by_aso_len:
            return score / len(aso)

        return score

    for flank in flank_windows:
        local_col = local_col_template.format(flank=flank)
        out_col = out_col_template.format(flank=flank)

        def compute_row(row):
            return compute_chimera(
                region=row.get(local_col, ""),
                aso=row.get(sequence_col, ""),
            )

        df[out_col] = df.apply(compute_row, axis=1)
        print(f"Finished calculating {out_col} for flank size {flank}")

    return df, suffix_array_cache

def add_global_chimera_feature_cds(
        df,
        sequence_col,
        cds_sequence_col="cds_sequence",
        out_col="chimera_score_global_CDS",
        step_size=1,
        normalize_by_aso_len=False,
        suffix_array_cache=None,
):
    """
    Adds a single global CDS chimera score:
      cds_sequence -> chimera_score_global_CDS

    Returns (df, suffix_array_cache).
    """
    if suffix_array_cache is None:
        suffix_array_cache = {}

    def safe_build_suffix_array(region_clean: str):
        if region_clean in suffix_array_cache:
            return suffix_array_cache[region_clean]
        sa = build_suffix_array(region_clean)
        suffix_array_cache[region_clean] = sa
        return sa

    def compute_row(row):
        region = row.get(cds_sequence_col, "")
        aso = row.get(sequence_col, "")

        if not isinstance(region, str) or not region.strip():
            return np.nan
        if not isinstance(aso, str) or not aso.strip():
            return np.nan

        region_clean = clean_region_from_aso(region, aso)
        if not isinstance(region_clean, str) or not region_clean.strip():
            return np.nan

        sa = safe_build_suffix_array(region_clean)
        score = calculate_chimera_ars(sa, aso, step_size=step_size)

        if normalize_by_aso_len:
            return score / len(aso)

        return score

    df[out_col] = df.apply(compute_row, axis=1)
    print(f"Finished calculating {out_col}")

    return df, suffix_array_cache
