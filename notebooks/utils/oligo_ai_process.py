import ast


def process_chemistry(mod_str):
    if not isinstance(mod_str, str) or not mod_str.startswith('['):
        return None, None

    try:
        mod_list = ast.literal_eval(mod_str)
        mod_set = {m.upper() for m in mod_list}

        # 1. Generate the CHEMICAL_PATTERN (e.g., MMMddd)
        pattern_chars = []
        for m in mod_list:
            m_up = m.upper()
            if m_up == 'MOE': pattern_chars.append('M')
            elif m_up == 'CET': pattern_chars.append('C')
            elif m_up == 'DNA': pattern_chars.append('d')
        pattern = "".join(pattern_chars)

        # 2. Determine the MODIFICATION label
        has_moe = 'MOE' in mod_set
        has_cet = 'CET' in mod_set or 'CET' in mod_set # handles case variations

        if has_moe and has_cet:
            label = "mixmer"
        elif has_moe:
            label = "MOE/5-methylcytosines/deoxy"
        elif has_cet:
            label = "cEt/5-methylcytosines/deoxy"
        else:
            label = "DNA"

        return pattern, label
    except:
        return None, None