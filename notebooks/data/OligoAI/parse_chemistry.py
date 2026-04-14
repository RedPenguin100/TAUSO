import ast

import pandas as pd

from tauso.data.consts import BACKBONE_MODS, CHEMICAL_PATTERN, MODIFICATION, PS_PATTERN, SUGAR_MODS


def _process_chemistry(mod_str):
    if not isinstance(mod_str, str) or not mod_str.startswith("["):
        return None, None

    try:
        mod_list = ast.literal_eval(mod_str)
        mod_set = {m.upper() for m in mod_list}

        # 1. Generate the CHEMICAL_PATTERN (e.g., MMMddd)
        pattern_chars = []
        for m in mod_list:
            m_up = m.upper()
            if m_up == "MOE":
                pattern_chars.append("M")
            elif m_up == "CET":
                pattern_chars.append("C")
            elif m_up == "DNA":
                pattern_chars.append("d")
        pattern = "".join(pattern_chars)

        # 2. Determine the MODIFICATION label
        has_moe = "MOE" in mod_set
        has_cet = "CET" in mod_set

        if has_moe and has_cet:
            label = "mixmer"
        elif has_moe:
            label = "MOE/5-methylcytosines/deoxy"
        elif has_cet:
            label = "cEt/5-methylcytosines/deoxy"
        else:
            label = "DNA"

        return pattern, label
    except (ValueError, SyntaxError):
        return None, None


def _parse_backbone(backbone_str):
    if not isinstance(backbone_str, str) or not backbone_str.startswith("["):
        return None

    try:
        raw_list = ast.literal_eval(backbone_str)
        mapping = {"PS": "*", "PO": "d", "<PAD>": ""}
        return "".join(mapping.get(item, "") for item in raw_list)
    except (ValueError, SyntaxError):
        return None


def assign_sugar(df: pd.DataFrame) -> pd.DataFrame:
    """Parses sugar chemistry strings and populates chemical pattern and modification labels."""
    if SUGAR_MODS not in df.columns:
        raise ValueError(f"column {SUGAR_MODS} is missing from the data, please populate it")

    # Assuming _process_chemistry is defined locally or imported
    results = df[SUGAR_MODS].apply(_process_chemistry)
    df[[CHEMICAL_PATTERN, MODIFICATION]] = pd.DataFrame(results.tolist(), index=df.index)
    return df


def assign_backbone(df: pd.DataFrame) -> pd.DataFrame:
    """Parses backbone chemistry strings and populates PS pattern labels."""
    if BACKBONE_MODS not in df.columns:
        raise ValueError(f"column {BACKBONE_MODS} is missing from the data, please populate it")

    # Assuming _parse_backbone is defined locally or imported
    df[PS_PATTERN] = df[BACKBONE_MODS].apply(_parse_backbone)
    return df


def assign_chemistry(df: pd.DataFrame) -> pd.DataFrame:
    """
    Parses sugar and backbone chemistry strings in a dataframe and populates new columns
    with their respective patterns and modification labels.
    """
    out_df = df.copy()

    out_df = assign_sugar(out_df)
    out_df = assign_backbone(out_df)

    return out_df
