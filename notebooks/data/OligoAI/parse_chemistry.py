import ast
import pandas as pd

from tauso.data.consts import CHEMICAL_PATTERN, MODIFICATION, PS_PATTERN


# --- Helper Functions ---
def _process_chemistry(mod_str):
    if not isinstance(mod_str, str) or not mod_str.startswith('['):
        return None, None

    try:
        mod_list = ast.literal_eval(mod_str)
        mod_set = {m.upper() for m in mod_list}

        # 1. Generate the CHEMICAL_PATTERN (e.g., MMMddd)
        pattern_chars = []
        for m in mod_list:
            m_up = m.upper()
            if m_up == 'MOE':
                pattern_chars.append('M')
            elif m_up == 'CET':
                pattern_chars.append('C')
            elif m_up == 'DNA':
                pattern_chars.append('d')
        pattern = "".join(pattern_chars)

        # 2. Determine the MODIFICATION label
        has_moe = 'MOE' in mod_set
        has_cet = 'CET' in mod_set

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
    if not isinstance(backbone_str, str) or not backbone_str.startswith('['):
        return None

    try:
        raw_list = ast.literal_eval(backbone_str)
        mapping = {
            'PS': '*',
            'PO': 'd',
            '<PAD>': ''
        }
        return "".join(mapping.get(item, '') for item in raw_list)
    except (ValueError, SyntaxError):
        return None


def populate_chemistry(df, sugar_col='sugar_mods', backbone_col='backbone_mods',
                       pattern_col=CHEMICAL_PATTERN, mod_col=MODIFICATION):
    """
    Parses sugar and backbone chemistry strings in a dataframe and populates new columns
    with their respective patterns and modification labels.
    """

    # --- DataFrame Application ---
    # Create a copy to avoid SettingWithCopyWarning if a slice was passed
    out_df = df.copy()

    # 1. Apply sugar modification logic
    if sugar_col in out_df.columns:
        results = out_df[sugar_col].apply(_process_chemistry)
        out_df[[pattern_col, mod_col]] = pd.DataFrame(results.tolist(), index=out_df.index)

    # 2. Apply backbone logic (if the column exists in your dataframe)
    if backbone_col in out_df.columns:
        out_df[PS_PATTERN] = out_df[backbone_col].apply(_parse_backbone)

    return out_df
