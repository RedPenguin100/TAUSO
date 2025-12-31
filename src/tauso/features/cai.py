
from __future__ import annotations
import re
import math
from typing import List, Union, Dict
from pathlib import Path
from typing import Iterable, List, Optional
import pandas as pd


STOP_CODONS = {"TAA", "TAG", "TGA"}

def _clean_dna(seq: str) -> str:
    """
    Normalize a nucleotide sequence:
      - Uppercase letters,
      - Convert RNA 'U' to DNA 'T',
      - Replace any non-ACGT character with 'N'.
    """
    seq = seq.upper().replace("U", "T")
    return re.sub(r"[^ACGT]", "N", seq)

# Keep the original family structure for backward compatibility and transparency.
FAMILIES: List[Dict[str, int]] = [
    {"TTT": 0, "TTC": 0},  # Phe
    {"ATT": 0, "ATC": 0, "ATA": 0},  # Ile
    {"GTT": 0, "GTC": 0, "GTA": 0, "GTG": 0},  # Val
    {"CCT": 0, "CCC": 0, "CCA": 0, "CCG": 0},  # Pro
    {"ACT": 0, "ACC": 0, "ACA": 0, "ACG": 0},  # Thr
    {"GCT": 0, "GCC": 0, "GCA": 0, "GCG": 0},  # Ala
    {"TAT": 0, "TAC": 0},  # Tyr
    {"CAT": 0, "CAC": 0},  # His
    {"CAA": 0, "CAG": 0},  # Gln
    {"AAT": 0, "AAC": 0},  # Asn
    {"AAA": 0, "AAG": 0},  # Lys
    {"GAT": 0, "GAC": 0},  # Asp
    {"GAA": 0, "GAG": 0},  # Glu
    {"TGT": 0, "TGC": 0},  # Cys
    {"GGT": 0, "GGC": 0, "GGA": 0, "GGG": 0},  # Gly
    {"TTA": 0, "TTG": 0, "CTT": 0, "CTC": 0, "CTA": 0, "CTG": 0},  # Leu
    {"TCT": 0, "TCC": 0, "TCA": 0, "TCG": 0, "AGT": 0, "AGC": 0},  # Ser
    {"CGT": 0, "CGC": 0, "CGA": 0, "CGG": 0, "AGA": 0, "AGG": 0},  # Arg
    {"ATG": 0},  # Met
    {"TGG": 0},  # Trp
]

def calc_CAI_weight(reference_seqs: Union[str, List[str]]):
    """
    Build CAI weights from one or more reference CDS sequences.

    Parameters
    ----------
    reference_seqs : str | list[str]
        One CDS string or a list of CDS strings (DNA alphabet: A/C/G/T).
        Typical usage: highly-expressed genes in the organism of interest.

    Returns
    -------
    weights_list : list[dict[str,float]]
        A list of dictionaries (one per synonymous family) where each codon
        weight is normalized by the family's maximum frequency.
        (Backward-compatible shape with your original code.)
    weights_flat : dict[str,float]
        A flat dictionary {codon: weight} for O(1) lookups.
    """
    if isinstance(reference_seqs, str):
        reference_seqs = [reference_seqs]

    # Fresh copy of families (counts start at zero)
    families: List[Dict[str, int]] = [{k: 0 for k in fam} for fam in FAMILIES]

    # Accumulate codon counts across all reference sequences
    for ref in reference_seqs:
        ref = _clean_dna(ref)
        for i in range(0, len(ref) - 2, 3):
            codon = ref[i:i+3]
            if "N" in codon or codon in STOP_CODONS:
                continue
            for fam in families:
                if codon in fam:
                    fam[codon] += 1
                    break

    # Normalize each family by its maximum count
    weights_list: List[Dict[str, float]] = []
    for fam in families:
        m = max(fam.values()) if fam else 0
        if m > 0:
            weights_list.append({k: (v / m) for k, v in fam.items()})
        else:
            # No information for this family; leave zeros (handled by epsilon later)
            weights_list.append({k: 0.0 for k in fam})

    # Build flat lookup
    weights_flat: Dict[str, float] = {}
    for fam in weights_list:
        weights_flat.update(fam)

    return weights_list, weights_flat

def calc_CAI(seq: str,
             weights: Union[List[Dict[str, float]], Dict[str, float]],
             epsilon: float = 1e-8) -> float:
    """
    Compute the Codon Adaptation Index (CAI) using a log-space geometric mean.

    Parameters
    ----------
    seq : str
        CDS sequence (DNA alphabet: A/C/G/T). Assumed to be in-frame;
        trailing bases are ignored. Internal stop codons are skipped.
    weights : list[dict[str,float]] | dict[str,float]
        Either the original "list of families" structure or a flat dict {codon: weight}.
    epsilon : float
        Small positive fallback when a codon weight is missing or zero.

    Returns
    -------
    float
        CAI value in [0, 1] (practically), or NaN if no valid codons remain.
    """
    clean = _clean_dna(seq)
    codons = [clean[i:i+3] for i in range(0, len(clean) - 2, 3)]
    # Ignore ambiguous codons and stops (same behavior as your original)
    codons = [c for c in codons if "N" not in c and c not in STOP_CODONS]
    if not codons:
        return float("nan")

    # Getter that works for either weights shape
    if isinstance(weights, list):
        def get_w(codon: str) -> float:
            for fam in weights:
                if codon in fam:
                    return fam[codon]
            return None
    else:
        def get_w(codon: str) -> float:
            return weights.get(codon)

    # Log-space geometric mean: exp( (1/N) * sum(log w_i) )
    s_log = 0.0
    N = 0
    for c in codons:
        w = get_w(c)
        if not w or w <= 0.0:
            w = epsilon
        s_log += math.log(w)
        N += 1

    return math.exp(s_log / max(N, 1))


############################################################################################################
# CAI reference weights from cell-line transcriptomes
DEFAULT_TRANSCRIPTOME_FILENAMES: List[str] = [
    "ACH-000232_transcriptome.csv",
    "ACH-000463_transcriptome.csv",
    "ACH-000681_transcriptome.csv",
    "ACH-000739_transcriptome.csv",
    "ACH-001086_transcriptome.csv",
    "ACH-001188_transcriptome.csv",
    "ACH-001328_transcriptome.csv",
]




def build_cai_reference_weights_from_transcriptomes(
    *,
    transcriptome_dir: Optional[Path] = None,
    transcriptome_filenames: Iterable[str] = DEFAULT_TRANSCRIPTOME_FILENAMES,
    top_n: int = 300,
    expr_col: str = "expression_norm",
    mutated_seq_col: str = "Mutated Transcript sequence",
    original_seq_col: str = "Original Transcript sequence",
) -> dict:
    """
    Build CAI reference weights from top-N transcript sequences by expression.

    Returns:
      weights_flat: dict[str, float]
        Codon -> CAI weight mapping used as the reference set.
    """
    if top_n <= 0:
        raise ValueError(f"top_n must be positive. Got: {top_n}")

    # Resolve default data path relative to this file (repo-stable)
    if transcriptome_dir is None:
        project_root = Path(__file__).resolve().parents[3]
        transcriptome_dir = project_root / "notebooks" / "transcripts" / "cell_line_expression"

    transcriptome_dir = Path(transcriptome_dir)

    files = [transcriptome_dir / fn for fn in transcriptome_filenames]
    missing = [p.name for p in files if not p.exists()]
    if missing:
        raise FileNotFoundError(f"Missing transcriptome files in {transcriptome_dir}: {missing}")

    dfs = [pd.read_csv(p) for p in files]
    transcript_df = pd.concat(dfs, ignore_index=True).drop_duplicates().copy()

    # Required columns
    for col in (expr_col, mutated_seq_col, original_seq_col):
        if col not in transcript_df.columns:
            raise KeyError(f"Missing required column '{col}' in transcriptome data.")

    # Prefer mutated, fallback to original
    transcript_df["ref_sequence"] = transcript_df[mutated_seq_col].fillna(
        transcript_df[original_seq_col]
    )

    # Top-N by expression
    ref_df = transcript_df.sort_values(expr_col, ascending=False).head(top_n)

    reference_seqs = ref_df["ref_sequence"].dropna().astype(str).tolist()
    if not reference_seqs:
        raise ValueError("No reference sequences found for CAI reference set.")

    # Build CAI weights
    _, weights_flat = calc_CAI_weight(reference_seqs)

    return weights_flat
