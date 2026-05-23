import numpy as np
import pandas as pd
from codonbias.scores import TrnaAdaptationIndex

# Human tRNA gene copy numbers, snapshot of GtRNAdb's Hsapi38 (hg38) summary
# fetched via codonbias.scores.fetch_GCN_from_GtRNAdb(genome="Hsapi38",
# domain="eukaryota"). Hardcoded so we do not depend on the network or on
# lxml at runtime. Refresh by re-running that helper if the source updates.
_HUMAN_TGCN = {
    "AAC": 9, "AAG": 9, "AAT": 15, "ACG": 7, "AGA": 9,
    "AGC": 26, "AGG": 9, "AGT": 9, "CAA": 6, "CAC": 13,
    "CAG": 9, "CAT": 20, "CCA": 7, "CCC": 5, "CCG": 4,
    "CCT": 5, "CGA": 4, "CGC": 4, "CGG": 4, "CGT": 5,
    "CTC": 8, "CTG": 13, "CTT": 15, "GAA": 10, "GAT": 3,
    "GCA": 29, "GCC": 14, "GCT": 8, "GTA": 13, "GTC": 13,
    "GTG": 9, "GTT": 25, "TAA": 4, "TAC": 5, "TAG": 3,
    "TAT": 5, "TCA": 1, "TCC": 9, "TCG": 6, "TCT": 6,
    "TGA": 4, "TGC": 8, "TGG": 7, "TGT": 6, "TTC": 8,
    "TTG": 6, "TTT": 12,
}

_TGCN_DF = pd.DataFrame(
    {"anti_codon": list(_HUMAN_TGCN.keys()), "GCN": list(_HUMAN_TGCN.values())}
)

_TAI_SCORER = TrnaAdaptationIndex(tGCN=_TGCN_DF, s_values="dosReis", genetic_code=1)


def compute_tAI(seq: str) -> float:
    if not isinstance(seq, str) or not seq.strip():
        return np.nan
    score = _TAI_SCORER.get_score(seq)
    return float(score) if np.isfinite(score) else np.nan
