from tauso.features.RNaseH_features import rnaseh1_dict, compute_rnaseh1_score


def add_RNaseH1_Krel(df, exp='R7_krel'):
    best_window_start_krel = {
        'R4a_krel': {10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 3, 18: 2, 19: 4, 20: 3, 21: 0, 22: 0, 25: 0},
        'R4b_krel': {10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 1, 18: 3, 19: 1, 20: 3, 21: 0, 22: 0, 25: 0},
        'R7_krel': {10: 0, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 3, 17: 2, 18: 4, 19: 6, 20: 4, 21: 0, 22: 0, 25: 0},
    }
    weights = rnaseh1_dict(exp)

    def score_row(row):
        length = len(row['Sequence'])
        pos = best_window_start_krel.get(exp, {}).get(length, 0)
        return compute_rnaseh1_score(row['Sequence'], weights, window_start=pos)

    col_name = f"RNaseH1_Krel_score_{exp}"
    df[col_name] = df.apply(score_row, axis=1)
    return df