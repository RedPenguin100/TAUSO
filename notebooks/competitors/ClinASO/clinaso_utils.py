"""
ClinASO competitor feature: the RNase H1 cleavage-preference PWM score.

Transcribed from ClinASO's `asodesigner/scr/_4.py` (vendored as the `ClinASO` submodule;
github.com/RedPenguin100/ClinASO, a fork of the Yunnan University platform published at
www.gapmerasodesign.com), using ClinASO's own position weight matrix
`asodesigner/text/human_rnaseH_pwm.txt` (4 nucleotides x 13 positions).
Chen et al., Mol Ther Nucleic Acids 37(2):102933, 2026 (doi:10.1016/j.omtn.2026.102933).

That PWM is byte-identical to `tauso.features.rnase_motifs.weights.R7`: both are the human
RNase H1 logFC weights from the R7 experiment of Kielpinski et al. 2017 (NAR; PMC5728404).
Negative logFC is the preferred (faster-cleaved) base there, which is why lower is better.
TAUSO's own `rnase_score_*` features apply the same weights but constrain the window to the
DNA gap; ClinASO instead sums four fixed offsets and ignores the gapmer design.

Algorithm, per `_4.py`: for offset in 0..3, take the 13-nt window ASO[2 + offset : 15 + offset],
reverse-complement it, and sum the PWM weights position-wise; the four windows sum to `total_pf`.
`_4.py` writes `total_pf` as the column ClinASO's downstream `_73.py` labels "RNase H score".

Sign convention: LOWER is better. ClinASO's final selection step `_8.py` sorts by
"RNase H score" ascending under its "efficiency" priority.

Two ways this score is not the whole of ClinASO. ClinASO is a rule-based designer, not a
trained efficacy model, and this PWM is its only potency-relevant per-ASO signal. And `_4.py`
discards candidates containing TCCA, GCTC, CCTG or GGG before scoring; here every ASO is
scored, so the column also covers candidates ClinASO would have rejected outright.

Windows span ASO[2:18], so ASOs shorter than 18 nt truncate the last windows. Because the
window is reverse-complemented before it is indexed against the PWM, a truncated window also
shifts the PWM register. Both behaviours match upstream.
"""
import functools
from pathlib import Path

import pandas as pd

from tauso.data.consts import ASO_SEQUENCE
from tauso.util import rna_to_dna

PWM_LENGTH = 13
_PWM_PATH = Path(__file__).resolve().parent / "ClinASO" / "asodesigner" / "text" / "human_rnaseH_pwm.txt"
_COMP_TABLE = str.maketrans("ATGCatgc", "TACGtacg")
_DEFAULT_SCORES = (0.0,) * PWM_LENGTH


def _dna_complement1(seq: str) -> str:
    """Reverse-complement, as ClinASO _4.py DNA_complement1."""
    return seq.translate(_COMP_TABLE)[::-1]


@functools.cache
def _load_pwm(path: Path = _PWM_PATH) -> dict[str, tuple[float, ...]]:
    dic: dict[str, tuple[float, ...]] = {}
    for line in Path(path).read_text().splitlines():
        parts = line.strip().split("\t")
        if len(parts) > 1:
            dic[parts[0]] = tuple(float(x) for x in parts[1 : 1 + PWM_LENGTH])
    if not dic:
        raise RuntimeError(
            f"ClinASO PWM empty/missing at {path}. "
            "Initialize the submodule: `git submodule update --init notebooks/competitors/ClinASO/ClinASO`."
        )
    return dic


def clinaso_total_pf(aso: str) -> float:
    """ClinASO RNase H1 preference score for a single ASO. Lower is better."""
    if not isinstance(aso, str) or not aso:
        return float("nan")
    dic_pwm = _load_pwm()
    sequence = rna_to_dna(aso)
    total_pf = 0.0
    for offset in range(4):
        start = 2 + offset
        sw = _dna_complement1(sequence[start : start + PWM_LENGTH])
        for j in range(min(len(sw), PWM_LENGTH)):
            total_pf += dic_pwm.get(sw[j], _DEFAULT_SCORES)[j]
    return total_pf


def populate_clinaso(data: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """Add the ClinASO RNase H1 PWM score column. Returns (df, feature_cols)."""
    df = data.copy()
    df["ClinASO_score"] = df[ASO_SEQUENCE].map(clinaso_total_pf)
    return df, ["ClinASO_score"]
