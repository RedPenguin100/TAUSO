# ================= RNA-mode rbp_compatibility =================
# Works entirely in RNA alphabet (A/C/G/U). No U->T conversions.
import numpy as np


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
