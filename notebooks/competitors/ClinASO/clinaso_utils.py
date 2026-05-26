"""
ClinASO competitor feature: the RNase H1 cleavage-preference PWM score ("total_pf").

Faithful transcription of ClinASO's scoring in notebooks/competitors/ClinASO/ClinASO/asodesigner/scr/_4.py
(github.com/wusuw/ClinASO; Chen et al., Mol Ther Nucleic Acids 37(2):102933, 2026), using
ClinASO's own position weight matrix (notebooks/competitors/ClinASO/ClinASO/asodesigner/text/human_rnaseH_pwm.txt,
4 nucleotides x 13 positions). ClinASO is a rule-based designer, not a trained efficacy model;
this PWM "total_pf" is its only potency-relevant per-ASO signal.

Algorithm (per _4.py): for offset in 0..3, take the 13-nt window ASO[2+offset : 2+offset+13],
reverse-complement it, sum the PWM weights position-wise; sum the 4 windows -> total_pf.

Follows the competitor populate pattern (cf. pfred_utils.populate_pfred / sfold_utils):
`populate_clinaso(data) -> (df, feature_cols)`.
"""
from pathlib import Path

import pandas as pd

from tauso.data.consts import SEQUENCE

PWM_LENGTH = 13
_PWM_PATH = Path(__file__).resolve().parent / "ClinASO" / "asodesigner" / "text" / "human_rnaseH_pwm.txt"
_COMP_TABLE = str.maketrans("ATGCatgc", "TACGtacg")
_DEFAULT_SCORES = (0.0,) * PWM_LENGTH
_DIC_PWM: dict | None = None


def _dna_complement1(seq: str) -> str:
    """Reverse-complement, identical to ClinASO _4.py DNA_complement1."""
    return seq.translate(_COMP_TABLE)[::-1]


def _load_pwm(path: Path = _PWM_PATH) -> dict:
    dic: dict = {}
    text = Path(path).read_text()
    for line in text.splitlines():
        parts = line.strip().split("\t")
        if len(parts) > 1:
            dic[parts[0]] = tuple(float(x) for x in parts[1:1 + PWM_LENGTH])
    if not dic:
        raise RuntimeError(
            f"ClinASO PWM empty/missing at {path}. "
            "Initialize the submodule: `git submodule update --init notebooks/competitors/ClinASO/ClinASO`."
        )
    return dic


def clinaso_total_pf(aso: str) -> float:
    """ClinASO RNase H1 preference score for a single ASO (faithful to _4.py)."""
    global _DIC_PWM
    if _DIC_PWM is None:
        _DIC_PWM = _load_pwm()
    if not isinstance(aso, str) or not aso:
        return float("nan")
    sequence = aso.upper().replace("U", "T")
    total_pf = 0.0
    for offset in range(4):
        start = 2 + offset
        sw = _dna_complement1(sequence[start:start + PWM_LENGTH])
        for j in range(min(len(sw), PWM_LENGTH)):
            total_pf += _DIC_PWM.get(sw[j], _DEFAULT_SCORES)[j]
    return total_pf


def populate_clinaso(data: pd.DataFrame):
    """
    Add the ClinASO RNase H1 PWM score column. Mirrors populate_pfred:
    returns (df, feature_cols).
    """
    df = data.copy()
    df["ClinASO_score"] = df[SEQUENCE].map(clinaso_total_pf)
    return df, ["ClinASO_score"]
