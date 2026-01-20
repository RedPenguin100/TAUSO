# ================= RNA-mode RBP_compatibility =================
# Works entirely in RNA alphabet (A/C/G/U). No U->T conversions.

import pandas as pd
from typing import List, Tuple, Dict, Optional
import pandas as pd
import math

def RBP_compatibility(aso_seq: str, target_rna_seq: str, rbp_motifs: List[str]) -> Optional[float]:
    """
    ASO–RBP compatibility (RNA alphabet):
      1.0 = fully compatible (no RBP clash)
      0.0 = strong clash with RBP motif
      None = if target_rna_seq empty

    Notes
    -----
    - Uses RNA complement map (A<->U, C<->G).
    - Locates the ASO binding site by searching the reverse-complement
      of the ASO (the target-side "sense" substring) within target_rna_seq.
      If not found, falls back to centering the ASO in the window.
    - For each exact-match motif occurrence, applies penalty = overlap / aso_len;
      final score = min over all overlaps.
    """

    if not target_rna_seq or len(str(target_rna_seq).strip()) == 0:
        return None

    def _norm_rna(s: str) -> str:
        # Uppercase, keep only A/C/G/U, strip spaces/newlines/tabs
        s = str(s).upper().replace(" ", "").replace("\t", "").replace("\n", "")
        # Remove anything not A/C/G/U
        return "".join(ch for ch in s if ch in {"A", "C", "G", "U"})

    aso = _norm_rna(aso_seq)
    target = _norm_rna(target_rna_seq)
    aso_len = len(aso)
    if aso_len == 0:
        return None

    # RNA reverse-complement
    comp = str.maketrans("ACGU", "UGCA")
    sense = aso.translate(comp)[::-1]

    # Locate ASO binding within the window; fallback = centered
    pos = target.find(sense)
    aso_start = pos if pos != -1 else max(0, (len(target) - aso_len) // 2)
    aso_end = aso_start + aso_len

    best_score = 1.0

    # Scan motifs (exact match) and compute maximal clash penalty
    for motif in rbp_motifs:
        m = _norm_rna(motif)
        if len(m) == 0 or len(m) > len(target):
            continue

        limit = len(target) - len(m) + 1
        for i in range(limit):
            if target[i:i+len(m)] == m:
                # overlap between ASO [aso_start, aso_end) and motif [i, i+len(m))
                overlap_start = max(aso_start, i)
                overlap_end   = min(aso_end,   i + len(m))
                overlap = max(0, overlap_end - overlap_start)
                if overlap > 0:
                    penalty = overlap / aso_len
                    best_score = min(best_score, 1 - penalty)

    return best_score





# ------------------------------------------------------------
# 2) Motif table loaders & cleaners (RNA mode)
# ------------------------------------------------------------

def _pick_col(df: pd.DataFrame, candidates: List[str]) -> str:
    """Pick a column by exact/lower match or substring fallback."""
    low = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in low:
            return low[cand.lower()]
    for c in df.columns:
        if any(cand.lower() in c.lower() for cand in candidates):
            return c
    raise KeyError(f"Could not find any of {candidates} in columns: {list(df.columns)}")

def _normalize_motif_series_rna(s: pd.Series) -> pd.Series:
    """Uppercase and keep only A/C/G/U."""
    return s.astype(str).str.upper().str.replace(r"[^ACGU]", "", regex=True)

def load_rbp_motif_table_rna(csv_path: str) -> pd.DataFrame:
    """
    Load motif CSV and normalize to RNA (A/C/G/U).
    Stores resolved column names in df.attrs['cols'] = {cell, motif, expr, role}.
    Expected columns (case-insensitive / fuzzy):
      - cell line / cell_line / cellline / cell / line
      - motif / sequence / motif_seq
      - expression / expr / tpm / fpkm / abundance
      - functional_role / role / category / class / group
    """
    df = pd.read_csv(csv_path)
    cell_col = _pick_col(df, ["cell line","cell_line","cellline","cell","line"])
    motif_col = _pick_col(df, ["motif","sequence","motif_seq"])
    expr_col  = _pick_col(df, ["expression","expr","tpm","fpkm","abundance"])
    role_col  = _pick_col(df, ["functional_role","role","category","class","group"])

    df = df.dropna(subset=[cell_col, motif_col, expr_col, role_col]).copy()
    df[motif_col] = _normalize_motif_series_rna(df[motif_col])

    df.attrs["cols"] = dict(cell=cell_col, motif=motif_col, expr=expr_col, role=role_col)
    return df


# ------------------------------------------------------------
# 3) Motif retrieval helpers (RNA mode)
# ------------------------------------------------------------

def get_relevant_motifs_rna(
    df: pd.DataFrame,
    *,
    cell_line: str,
    roles: List[str],
    expr_threshold: float
) -> List[str]:
    """
    Return unique RNA motifs (A/C/G/U) for a given `cell_line` and list of `roles`,
    filtered by Expression >= expr_threshold.
    `df` must be the result of `load_rbp_motif_table_rna(...)`.
    """
    cols = df.attrs["cols"]
    sub = df[
        (df[cols["cell"]].astype(str) == str(cell_line)) &
        (df[cols["role"]].isin([str(r) for r in roles])) &
        (df[cols["expr"]].astype(float) >= float(expr_threshold))
    ]
    return sub[cols["motif"]].dropna().astype(str).unique().tolist()

def get_motifs_by_cellline_rna(
    df: pd.DataFrame,
    *,
    roles: List[str],
    expr_threshold: float
) -> Dict[str, List[str]]:
    """
    For every cell line present in `df`, return its relevant motif list
    (RNA A/C/G/U) given `roles` and `Expression >= expr_threshold`.
    """
    cols = df.attrs["cols"]
    out: Dict[str, List[str]] = {}
    for cl in df[cols["cell"]].astype(str).unique():
        out[cl] = get_relevant_motifs_rna(
            df, cell_line=cl, roles=roles, expr_threshold=expr_threshold
        )
    return out


###################################################################################
def overlap_density(aso_seq: str, target_seq: str, motifs: list[str]) -> float | None:
    """
    Calculate the RBP overlap density between an ASO and a target RNA sequence.

    The function measures what fraction of the ASO binding site is covered 
    by RBP motifs (for a given role, e.g., destabilizer). It does so by:
      1. Normalizing sequences to RNA alphabet (A/C/G/U).
      2. Locating the ASO binding site in the target (using reverse complement).
      3. Identifying all motif occurrences that overlap with the ASO site.
      4. Merging overlaps into a union to avoid double counting.
      5. Returning the fraction of the ASO length covered by these overlaps.

    Returns:
      - A float between 0.0 and 1.0:
          0.0 → no overlap (ASO not covered by motifs)
          1.0 → ASO fully covered by motifs
          intermediate values → partial coverage
      - None if either ASO or target sequence is empty.
    """

    # 1. Normalize sequences to RNA alphabet
    def norm(s: str) -> str:
        s = str(s).upper().replace(" ", "").replace("\t", "").replace("\n", "")
        return "".join(ch for ch in s if ch in {"A", "C", "G", "U"})
    aso = norm(aso_seq)
    target = norm(target_seq)
    if not target or not aso:
        return None

    # 2. Locate ASO binding site in the target
    comp = str.maketrans("ACGU", "UGCA")       # RNA complement mapping
    sense = aso.translate(comp)[::-1]          # reverse complement of the ASO
    pos = target.find(sense)                   # search within the target sequence
    if pos == -1:                              # fallback if not found
        pos = max(0, (len(target) - len(aso)) // 2)
    aso_start, aso_end = pos, pos + len(aso)

    # 3. Find motif occurrences that overlap with the ASO
    intervals = []
    for m in motifs:
        mm = norm(m)
        if not mm or len(mm) > len(target):
            continue
        limit = len(target) - len(mm) + 1
        for i in range(limit):
            if target[i:i+len(mm)] == mm:      # motif occurrence
                a = max(aso_start, i)          # overlap with ASO interval
                b = min(aso_end,   i + len(mm))
                if b > a:
                    intervals.append((a, b))

    if not intervals:
        return 0.0

    # 4. Merge intervals (union) to measure net coverage
    intervals.sort()
    total = 0
    cur_s, cur_e = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_e:                         # overlapping intervals
            cur_e = max(cur_e, e)
        else:                                  # close current and start new
            total += (cur_e - cur_s)
            cur_s, cur_e = s, e
    total += (cur_e - cur_s)

    # 5. Compute overlap density (fraction of ASO covered)
    return total / len(aso) if len(aso) > 0 else None
 
 ##############################################################################################################
def add_role_contrast_features(filtered, role_a="destabilizer", role_b="stabilizer", windows=None):
    """
    Add delta and ratio features comparing overlap densities of two roles.

    Parameters
    ----------
    filtered : pd.DataFrame
        DataFrame that already contains per-window and global overlap density columns.
    role_a : str
        Numerator role (e.g., "destabilizer").
    role_b : str
        Denominator role (e.g., "stabilizer").
    windows : list[int] or None
        List of window sizes (CDS_WINDOWS). If None, only global features are added.

    Returns
    -------
    list of str
        Names of the newly created columns.
    """

    # clean suffix for column names
    def clean(role: str) -> str:
        return role.replace("/", "_").replace(" ", "_")

    suf_a = clean(role_a)
    suf_b = clean(role_b)

    created_cols = []

    if windows is not None:
        for flank in windows:
            col_a = f"RBP_overlap_density_{role_a}_{flank}"
            col_b = f"RBP_overlap_density_{role_b}_{flank}"

            delta_col = f"RBP_overlap_delta_{suf_a}_vs_{suf_b}_{flank}"
            ratio_col = f"RBP_overlap_ratio_{suf_a}_vs_{suf_b}_{flank}"

            filtered[delta_col] = filtered[col_a] - filtered[col_b]
            filtered[ratio_col] = filtered[col_a] / (filtered[col_b] + 1e-6)

            created_cols.extend([delta_col, ratio_col])

    # Global
    col_a = f"RBP_overlap_density_{role_a}_global"
    col_b = f"RBP_overlap_density_{role_b}_global"

    delta_col = f"RBP_overlap_delta_{suf_a}_vs_{suf_b}_global"
    ratio_col = f"RBP_overlap_ratio_{suf_a}_vs_{suf_b}_global"

    filtered[delta_col] = filtered[col_a] - filtered[col_b]
    filtered[ratio_col] = filtered[col_a] / (filtered[col_b] + 1e-6)

    created_cols.extend([delta_col, ratio_col])

    return created_cols


##############################################################################################################
def positional_splits_exact_rna(aso_seq: str, target_seq: str, motifs: list[str]):
    """
    Compute (left, core, right) fractions of the ASO covered by RBP motifs
    using exact RNA matches.

    Steps:
    1. Normalize ASO and target sequences to RNA alphabet (A/C/G/U).
    2. Locate the ASO binding site on the target via RNA reverse-complement.
    3. Find all exact motif matches in the target and collect their overlaps
       with the ASO binding site.
    4. Merge overlapping intervals into disjoint segments.
    5. Split the ASO into three equal parts (left, core, right) and calculate
       what fraction of the ASO length is covered in each part.

    Parameters
    ----------
    aso_seq : str
        The ASO sequence (RNA, A/C/G/U).
    target_seq : str
        The target sequence (CDS or window context, RNA).
    motifs : list of str
        List of RBP motifs (RNA alphabet).

    Returns
    -------
    tuple of floats (left_frac, core_frac, right_frac)
        Fractions of ASO nucleotides covered by motif overlaps in each region.
        Values are in [0,1].
        If no overlaps exist, returns (0.0, 0.0, 0.0).
    """

    # --- normalize to RNA alphabet ---
    def norm(s: str) -> str:
        s = str(s).upper().replace(" ", "").replace("\t", "").replace("\n", "")
        return "".join(ch for ch in s if ch in {"A", "C", "G", "U"})

    aso = norm(aso_seq)
    target = norm(target_seq)
    if not aso or not target:
        return (0.0, 0.0, 0.0)

    # --- locate ASO on target (via reverse-complement in RNA) ---
    comp = str.maketrans("ACGU", "UGCA")
    sense = aso.translate(comp)[::-1]
    pos = target.find(sense)
    if pos == -1:
        # fallback: center the ASO if not found
        pos = max(0, (len(target) - len(aso)) // 2)
    aso_start, aso_end = pos, pos + len(aso)

    # --- collect overlap intervals with motifs ---
    intervals = []
    for m in motifs:
        mm = norm(m)
        if not mm or len(mm) > len(target):
            continue
        limit = len(target) - len(mm) + 1
        for i in range(limit):
            if target[i:i+len(mm)] == mm:
                a = max(aso_start, i)
                b = min(aso_end, i + len(mm))
                if b > a:
                    intervals.append((a, b))

    if not intervals:
        return (0.0, 0.0, 0.0)

    # --- merge intervals ---
    intervals.sort()
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))

    # --- split into thirds and compute fractions ---
    thirds = [0, len(aso) // 3, 2 * len(aso) // 3, len(aso)]

    def part_len(a, b, L, R):
        return max(0, min(b, R) - max(a, L))

    L1 = L2 = L3 = 0
    for s, e in merged:
        L1 += part_len(s, e, thirds[0], thirds[1])
        L2 += part_len(s, e, thirds[1], thirds[2])
        L3 += part_len(s, e, thirds[2], thirds[3])

    denom = max(1, len(aso))
    return L1 / denom, L2 / denom, L3 / denom


##############################################################################################################
def get_motifs_with_expr_by_cellline_rna(
    df: pd.DataFrame,
    *,
    roles: List[str],
    expr_threshold: float
) -> Dict[str, List[Tuple[str, float]]]:
    """
    For every cell line present in `df`, return a list of (motif, expression)
    filtered by roles and Expression >= expr_threshold.
    `df` must come from load_rbp_motif_table_rna(...).
    """
    cols = df.attrs["cols"]
    out: Dict[str, List[Tuple[str, float]]] = {}
    for cl in df[cols["cell"]].astype(str).unique():
        sub = df[
            (df[cols["cell"]].astype(str) == str(cl)) &
            (df[cols["role"]].isin([str(r) for r in roles])) &
            (pd.to_numeric(df[cols["expr"]], errors="coerce") >= float(expr_threshold))
        ][[cols["motif"], cols["expr"]]].dropna()
        # normalize motif to RNA alphabet and coerce expr to float
        motifs = []
        for _, row in sub.iterrows():
            m = str(row[cols["motif"]]).upper()
            m = "".join(ch for ch in m if ch in {"A","C","G","U"})
            if not m:
                continue
            try:
                e = float(row[cols["expr"]])
            except Exception:
                continue
            motifs.append((m, e))
        out[str(cl)] = motifs
    return out
##########################################################################################3
from typing import Tuple

def expression_weighted_density(
    aso_seq: str,
    target_seq: str,
    motifs_with_expr: List[Tuple[str, float]],
) -> Optional[float]:
    """
    Expression-weighted RBP overlap density.
    Sums overlap_length * expression over all motif overlaps, normalized by ASO length.
    Returns 0.0 if no overlaps; None if sequences empty.
    """

    def norm(s: str) -> str:
        s = str(s).upper().replace(" ", "").replace("\t", "").replace("\n", "")
        return "".join(ch for ch in s if ch in {"A","C","G","U"})

    aso = norm(aso_seq)
    target = norm(target_seq)
    if not aso or not target:
        return None

    comp = str.maketrans("ACGU", "UGCA")
    sense = aso.translate(comp)[::-1]
    pos = target.find(sense)
    if pos == -1:
        pos = max(0, (len(target) - len(aso)) // 2)
    aso_start, aso_end = pos, pos + len(aso)

    total = 0.0
    for motif, expr in motifs_with_expr:
        mm = norm(motif)
        if not mm or len(mm) > len(target):
            continue
        limit = len(target) - len(mm) + 1
        for i in range(limit):
            if target[i:i+len(mm)] == mm:
                a = max(aso_start, i)
                b = min(aso_end,   i + len(mm))
                if b > a:
                    total += (b - a) * float(expr)

    return total / len(aso) if len(aso) > 0 else None
####################################################################################
from typing import Dict, List, Tuple
import pandas as pd

def get_motif_identity_pairs_by_cellline(
    df: pd.DataFrame,
    role: str,
    expr_threshold: float,
    identity_col: str = "Gene_name"
) -> Dict[str, List[Tuple[str, str]]]:
    """
    Build mapping: cell_line -> list of (motif_RNA, identity).
    Filters by role and Expression >= expr_threshold.
    Uses 'Gene_name' (or chosen identity_col) as identity.
    """
    out: Dict[str, List[Tuple[str, str]]] = {}
    for cl in df["cell line"].astype(str).unique():
        sub = df.loc[
            (df["cell line"].astype(str) == cl) &
            (df["functional_role"] == role) &
            (pd.to_numeric(df["Expression"], errors="coerce") >= float(expr_threshold)),
            ["Motif", identity_col]
        ].dropna()

        pairs: List[Tuple[str, str]] = []
        for m, ident in sub.itertuples(index=False):
            motif_rna = "".join(ch for ch in str(m).upper() if ch in {"A","C","G","U"})
            if motif_rna:
                pairs.append((motif_rna, str(ident)))
        out[str(cl)] = pairs
    return out
####################################################################################
def motif_diversity_window_ratio(window_seq: str, motif_identity_pairs: list[tuple[str, str]]) -> float:
    """
    Compute the motif diversity ratio within a given sequence window.

    Definition:
        Diversity ratio = (# of unique Gene_name identities among all motif hits in the window)
                          / (total number of motif hits in the window)

    Inputs:
        window_seq : str
            Nucleotide sequence (the window around the ASO binding site, or the full CDS).
        motif_identity_pairs : list of (motif, identity)
            Each tuple contains a motif sequence (string of A/C/G/U) and the associated identity
            (typically Gene_name of the RBP that recognizes that motif).

    Procedure:
        - Normalize sequences to RNA alphabet (A,C,G,U).
        - Scan the entire window_seq for each motif.
        - For every occurrence of a motif, record its identity (e.g. the RBP name).
        - Count the total number of hits and the number of unique identities.

    Output:
        float
            Diversity ratio in [0,1].
            Returns 0.0 if no motif hits are found in the window.

    Notes:
        - If all hits come from a single RBP identity, ratio will be low (~0).
        - If each hit comes from a different RBP identity, ratio will approach 1.0.
        - This metric captures how many distinct RBPs are represented relative to
          the total number of motif matches in the context window.
    """
    def norm(s: str) -> str:
        s = str(s).upper().replace(" ", "").replace("\t", "").replace("\n", "")
        return "".join(ch for ch in s if ch in {"A","C","G","U"})

    win = norm(window_seq)
    if not win:
        return 0.0

    ids = []
    for motif, gene in motif_identity_pairs:
        mm = norm(motif)
        if not mm or len(mm) > len(win):
            continue
        for i in range(len(win) - len(mm) + 1):
            if win[i:i+len(mm)] == mm:
                ids.append(str(gene))

    if not ids:
        return 0.0
    return len(set(ids)) / len(ids)

####################################################################################
from typing import List, Tuple, Optional

def dominant_rbp_fraction(
    aso_seq: str,
    target_seq: str,
    motif_identity_pairs: List[Tuple[str, str]],
) -> Optional[float]:
    """
    Dominant RBP identity fraction over the ASO binding site.

    Interpretation:
      - 1.0  : all overlapping RBP coverage is dominated by a single RBP identity
      - 0.0  : no motif overlaps at all
      - None : missing/empty ASO or target sequence

    Parameters
    ----------
    aso_seq : str
        ASO sequence in RNA alphabet (A/C/G/U).
    target_seq : str
        Target RNA sequence (window or CDS) in RNA alphabet.
    motif_identity_pairs : list[(motif, identity)]
        motif: RNA motif sequence (A/C/G/U)
        identity: e.g., Gene_name of the RBP

    Returns
    -------
    float in [0,1] or None
    """

    def norm(s: str) -> str:
        s = str(s).upper().replace(" ", "").replace("\t", "").replace("\n", "")
        return "".join(ch for ch in s if ch in {"A", "C", "G", "U"})

    aso = norm(aso_seq)
    target = norm(target_seq)
    if not aso or not target:
        return None

    # locate ASO binding site in target using RNA reverse-complement
    comp = str.maketrans("ACGU", "UGCA")
    sense = aso.translate(comp)[::-1]
    pos = target.find(sense)
    if pos == -1:
        pos = max(0, (len(target) - len(aso)) // 2)

    aso_start, aso_end = pos, pos + len(aso)

    # accumulate overlap coverage per identity
    coverage_by_identity = {}

    for motif, identity in motif_identity_pairs:
        mm = norm(motif)
        if not mm or len(mm) > len(target):
            continue

        for i in range(len(target) - len(mm) + 1):
            if target[i:i + len(mm)] == mm:
                a = max(aso_start, i)
                b = min(aso_end, i + len(mm))
                if b > a:
                    coverage_by_identity[identity] = (
                        coverage_by_identity.get(identity, 0) + (b - a)
                    )

    if not coverage_by_identity:
        return 0.0

    total = sum(coverage_by_identity.values())
    dominant = max(coverage_by_identity.values())
    return

