from __future__ import annotations
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[3]  # TAUSO
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd

from src.tauso.genome.read_human_genome import get_locus_to_data_dict
from src.tauso.consts import CACHE_DIR
from pathlib import Path
from src.tauso.util import get_antisense as _default_get_antisense


PROJECT_ROOT = Path(__file__).resolve().parents[3]
DATA_MRNA_PATH = PROJECT_ROOT / "notebooks" / "transcripts" / "cell_line_expression"

assert DATA_MRNA_PATH.exists(), f"Expected transcriptome dir not found: {DATA_MRNA_PATH}"


# =============================================================================
# Defaults: update once here, then never touch scripts again
# =============================================================================

DEFAULT_TRANSCRIPTOME_FILENAMES: List[str] = [
    # Original transcriptome (full-length)
    "ACH-000232_transcriptome.csv",
    "ACH-000463_transcriptome.csv",
    "ACH-000681_transcriptome.csv",
    "ACH-000739_transcriptome.csv",
    "ACH-001086_transcriptome.csv",
    "ACH-001188_transcriptome.csv",
    "ACH-001328_transcriptome.csv",

]

DEFAULT_MRNA_RELATIVE_DIR = Path("notebooks") / "transcripts" / "cell_line_expression"


# =============================================================================
# Internal helpers: sequences
# =============================================================================

def _norm_rna_to_dna(seq: str) -> str:
    """Normalize RNA to DNA alphabet and strip whitespace."""
    return str(seq).upper().replace("U", "T").replace(" ", "").replace("\t", "").replace("\n", "")


def _to_str_seq(x) -> str:
    """
    Coerce sequence-like (list/np.array/Series) or string to a clean uppercase DNA string.
    Ensures slicing returns a plain string and avoids pandas iterable assignment bugs.
    """
    if isinstance(x, str):
        s = x
    else:
        try:
            s = "".join(list(x))
        except Exception:
            s = str(x)
    return s.replace(" ", "").replace("\t", "").replace("\n", "").replace("U", "T").upper()


def _build_spliced_mrna_from_exons(pre_mrna: str, exon_indices: Iterable[Tuple[int, int]]) -> str:
    """
    Build exon-joined mRNA by concatenating exon slices out of pre_mrna.

    Assumption (same as your current code):
      - pre_mrna corresponds to genomic strand and starts at exon_indices[0][0]
      - exon intervals are used directly
    """
    exon_indices = list(exon_indices)
    if not exon_indices:
        return ""
    pre_genome_start = exon_indices[0][0]
    parts = []
    for exon_start, exon_end in exon_indices:
        pm_start = exon_start - pre_genome_start
        pm_end = exon_end - pre_genome_start
        parts.append(pre_mrna[pm_start:pm_end])
    return "".join(parts)


# =============================================================================
# gene_to_data caching: this wraps your existing snippet
# =============================================================================

def load_or_build_gene_to_data(
    *,
    gene_subset: Optional[Iterable[str]] = None,
    include_introns: bool = True,
    cache_name: str = "gene_to_data_simple_cache.pickle",
) -> Dict[str, object]:
    """
    Inputs you may pass:
      - gene_subset: iterable of gene symbols (recommended to keep cache smaller and faster)
      - include_introns: True keeps introns in full_mrna
      - cache_name: cache filename under CACHE_DIR

    Output you get:
      - gene_to_data: dict[str, locus_info]
        locus_info is produced by get_locus_to_data_dict and must include:
          - full_mrna
          - exon_indices
          - cds_start, cds_end
    """
    cache_path = CACHE_DIR / cache_name
    if cache_path.exists():
        import pickle
        with open(cache_path, "rb") as f:
            return pickle.load(f)

    gene_to_data = get_locus_to_data_dict(include_introns=include_introns, gene_subset=gene_subset)
    import pickle
    with open(cache_path, "wb") as f:
        pickle.dump(gene_to_data, f)
    return gene_to_data


# =============================================================================
# External mRNA loader: manual file list, with partial-mutated awareness
# =============================================================================

def _load_mrna_by_gene_from_files(files: List[Path], seq_column: str) -> Dict[str, str]:
    """
    Load {Gene -> sequence} from CSV files.
    Keeps the longest clean A/C/G/T sequence per gene.
    """
    rows = []
    for f in files:
        df = pd.read_csv(f, usecols=["Gene", seq_column])
        df[seq_column] = df[seq_column].map(_norm_rna_to_dna)
        df = df[df[seq_column].str.fullmatch(r"[ACGT]+", na=False)]
        rows.append(df)

    if not rows:
        return {}

    big = pd.concat(rows, ignore_index=True)
    big["len"] = big[seq_column].str.len()
    chosen = big.sort_values(["Gene", "len"], ascending=[True, False]).drop_duplicates("Gene")
    return dict(zip(chosen["Gene"], chosen[seq_column]))


def _build_gene_to_mrna_real(
    data_path: Path,
    filenames: List[str],
    *,
    prefer_mutated: bool = True,
    allow_partial_mutated_override: bool = False,
) -> Tuple[Dict[str, str], Dict[str, dict]]:
    """
    Builds:
      - gene_to_mrna_real: {Gene -> chosen sequence}
      - gene_to_mrna_meta: {Gene -> metadata about choice}

    Important behavior:
      - "*_transcriptome.csv" is treated as full-length.
      - "*_mutated_transcriptome_premRNA.100.csv" is treated as partial (not full-length).

    Default safety:
      - allow_partial_mutated_override=False means partial mutated does NOT override full-length original.
      - If you set allow_partial_mutated_override=True, partial mutated can override.

    Returns:
      gene_to_mrna_real and gene_to_mrna_meta (source, is_full_length, availability flags).
    """
    files = [data_path / fn for fn in filenames]
    missing = [p.name for p in files if not p.exists()]
    if missing:
        raise FileNotFoundError(f"Missing transcriptome files in {data_path}: {missing}")

    mutated_re = re.compile(r"_mutated_transcriptome_premRNA\.100\.csv$", re.IGNORECASE)

    orig_files = [
        p for p in files
        if p.name.endswith("_transcriptome.csv") and (mutated_re.search(p.name) is None)
    ]
    mut_files = [p for p in files if mutated_re.search(p.name) is not None]

    d_orig: Dict[str, str] = {}
    if orig_files:
        try:
            d_orig = _load_mrna_by_gene_from_files(orig_files, "Original Transcript sequence")
        except ValueError:
            d_orig = _load_mrna_by_gene_from_files(orig_files, "Original Transcript Sequence")

    d_mut: Dict[str, str] = {}
    if mut_files:
        try:
            d_mut = _load_mrna_by_gene_from_files(mut_files, "Mutated Transcript sequence")
        except ValueError:
            d_mut = _load_mrna_by_gene_from_files(mut_files, "Mutated Transcript Sequence")

    gene_to_mrna_real: Dict[str, str] = {}
    gene_to_mrna_meta: Dict[str, dict] = {}

    all_genes = set(d_orig) | set(d_mut)
    for g in all_genes:
        orig_seq = d_orig.get(g)
        mut_seq = d_mut.get(g)

        chosen_seq = None
        source = None
        is_full_length = None

        if prefer_mutated and (mut_seq is not None) and allow_partial_mutated_override:
            chosen_seq = mut_seq
            source = "mutated_premRNA100"
            is_full_length = False
        elif orig_seq is not None:
            chosen_seq = orig_seq
            source = "original_transcriptome"
            is_full_length = True
        elif mut_seq is not None:
            chosen_seq = mut_seq
            source = "mutated_premRNA100"
            is_full_length = False

        if chosen_seq is None:
            continue

        gene_to_mrna_real[g] = chosen_seq
        gene_to_mrna_meta[g] = {
            "source": source,
            "is_full_length": is_full_length,
            "has_original": orig_seq is not None,
            "has_mutated_premRNA100": mut_seq is not None,
        }

    return gene_to_mrna_real, gene_to_mrna_meta


# =============================================================================
# Public API: one function you call from every script
# =============================================================================

def add_external_mrna_and_context_columns(
    df: pd.DataFrame,
    *,
    # Required: you must pass these
    sequence_col: str,
    canonical_gene_col: str,
    get_antisense_fn = _default_get_antisense,

    # Recommended: keeps runtime and cache smaller
    gene_subset: Optional[Iterable[str]] = None,

    # Behavior toggles
    flank_sizes_premrna: Iterable[int] = (20, 30, 40, 50, 60, 70),
    flank_sizes_cds: Iterable[int] = (20, 30, 40, 50, 60, 70),
    add_region_is_local_flags: bool = True,
    filter_not_found: bool = True,
    inplace: bool = False,

    # External mRNA selection behavior
    prefer_mutated: bool = True,
    allow_partial_mutated_override: bool = False,

    # External mRNA metadata columns (recommended for transparency)
    store_mrna_choice_columns: bool = True,

    # File discovery controls (manual list is default)
    transcriptome_filenames: List[str] = DEFAULT_TRANSCRIPTOME_FILENAMES,
    mrna_relative_dir: Path = DEFAULT_MRNA_RELATIVE_DIR,
    project_root: Optional[Path] = None,

    # gene_to_data cache control
    gene_to_data_cache_name: str = "gene_to_data_simple_cache.pickle",
) -> pd.DataFrame:
    """
    What you must provide:
      - df: your ASO DataFrame
      - sequence_col: the column name in df that contains the ASO antisense sequence
      - canonical_gene_col: the column name in df that contains gene symbols (must match gene_to_data keys)
      - get_antisense_fn: function that converts antisense ASO to its sense complement for searching pre_mRNA

    Defaults already handled inside:
      - gene_to_data is loaded from CACHE_DIR/gene_to_data_cache_name
        or built via get_locus_to_data_dict(include_introns=True, gene_subset=gene_subset)
      - transcriptome files are taken from transcriptome_filenames (manual list)
      - transcriptome base directory is <project_root>/scripts/data_genertion/cell_line_expression
        If project_root is None, it is inferred relative to this file path.

    What you get back:
      A DataFrame with these added columns:
        - sense_start: index of binding site in pre-mRNA (or -1 if not found)
        - sense_length: ASO length
        - sense_type: "exon" or "intron"
        - flank_sequence_{fs}: pre-mRNA window around site
        - cds_sequence: CDS string for the gene
        - in_coding_region: True if site maps into CDS
        - local_coding_region_around_ASO_{fs}: CDS window around site (only when in CDS)
        - region_is_local_{fs}: 0/1 flag whether local CDS window exists (optional)

      If store_mrna_choice_columns=True, also:
        - mrna_chosen_source: "original_transcriptome" or "mutated_premRNA100"
        - mrna_chosen_is_full_length: 1 for full-length, 0 for partial

    Important note:
      - mutated_transcriptome_premRNA.100 is partial. Default is NOT to let it override full-length original.
        If you really want partial mutated to override, set allow_partial_mutated_override=True.
    """
    if not inplace:
        df = df.copy()

    if gene_subset is not None:
        gene_subset_set = set(gene_subset)
        df = df[df[canonical_gene_col].isin(gene_subset_set)].copy()
    else:
        gene_subset_set = None
    """""
    # Load gene_to_data from cache (or build if missing)
    gene_to_data = load_or_build_gene_to_data(
        gene_subset=gene_subset_set,
        include_introns=True,
        cache_name=gene_to_data_cache_name,
    )
    """
    gene_to_data = get_locus_to_data_dict(
    include_introns=True,
    gene_subset=gene_subset_set
)

    # Resolve project root
    if project_root is None:
        # Adjust parents[k] once if needed, then keep stable.
        this_dir = Path(__file__).resolve().parent
        project_root = this_dir.parents[3]
    data_path = DATA_MRNA_PATH

    # Load external transcriptome sequences
    gene_to_mrna_real, gene_to_mrna_meta = _build_gene_to_mrna_real(
        data_path=data_path,
        filenames=transcriptome_filenames,
        prefer_mutated=prefer_mutated,
        allow_partial_mutated_override=allow_partial_mutated_override,
    )

    # Column names created
    SENSE_START = "sense_start"
    SENSE_LENGTH = "sense_length"
    SENSE_TYPE = "sense_type"
    CDS_SEQUENCE = "cds_sequence"
    IN_CODING_REGION = "in_coding_region"

    flank_sizes_premrna = list(flank_sizes_premrna)
    flank_sizes_cds = list(flank_sizes_cds)

    # Initialize columns
    df[SENSE_START] = -1
    df[SENSE_LENGTH] = 0
    df[SENSE_TYPE] = "NA"
    df[CDS_SEQUENCE] = ""
    df[IN_CODING_REGION] = False

    for fs in flank_sizes_premrna:
        df[f"flank_sequence_{fs}"] = ""
    for fs in flank_sizes_cds:
        df[f"local_coding_region_around_ASO_{fs}"] = ""

    if store_mrna_choice_columns:
        df["mrna_chosen_source"] = ""
        df["mrna_chosen_is_full_length"] = 0

    # Cache CDS per gene for speed: (cds_seq, genome_to_cds_map)
    gene_to_cds_info: Dict[str, Tuple[str, Dict[int, int]]] = {}

    for ridx, row in df.iterrows():
        gene_name = row[canonical_gene_col]
        locus_info = gene_to_data.get(gene_name)
        if locus_info is None:
            continue

        # Update per-row metadata about external mRNA choice (transparent, avoids self-deception)
        if store_mrna_choice_columns:
            meta = gene_to_mrna_meta.get(gene_name)
            if meta is not None:
                df.at[ridx, "mrna_chosen_source"] = meta["source"]
                df.at[ridx, "mrna_chosen_is_full_length"] = int(bool(meta["is_full_length"]))

        # Keep using your existing pre-mRNA for site finding and flanks
        pre_mrna = _to_str_seq(locus_info.full_mrna)

        antisense = _to_str_seq(row[sequence_col])
        sense = _to_str_seq(get_antisense_fn(antisense))

        # Locate site on pre-mRNA
        idx = pre_mrna.find(sense)
        df.at[ridx, SENSE_START] = idx
        df.at[ridx, SENSE_LENGTH] = len(antisense)

        if idx == -1:
            continue

        # Genomic correction (kept as-is)
        genome_corrected_index = idx + locus_info.exon_indices[0][0]

        # exon / intron classification (kept as-is)
        region_type = "intron"
        for exon_start, exon_end in locus_info.exon_indices:
            if exon_start <= genome_corrected_index <= exon_end:
                region_type = "exon"
                break
        df.at[ridx, SENSE_TYPE] = region_type

        # pre-mRNA flanks
        for fs in flank_sizes_premrna:
            flank_start = max(0, idx - fs)
            flank_end = min(len(pre_mrna), idx + len(sense) + fs)
            df.at[ridx, f"flank_sequence_{fs}"] = pre_mrna[flank_start:flank_end]

        # Build CDS + genome->CDS map (cached per gene)
        if gene_name not in gene_to_cds_info:
            cds_chars: List[str] = []
            genome_to_cds_map: Dict[int, int] = {}

            mrna_idx = 0
            for exon_start, exon_end in locus_info.exon_indices:
                for gpos in range(exon_start, exon_end):
                    if mrna_idx >= len(pre_mrna):
                        break
                    if locus_info.gene_start <= gpos <= locus_info.gene_end:
                        cds_chars.append(pre_mrna[mrna_idx])
                        genome_to_cds_map[gpos] = len(cds_chars) - 1
                    mrna_idx += 1

            cds_seq = "".join(cds_chars)
            gene_to_cds_info[gene_name] = (cds_seq, genome_to_cds_map)

        cds_seq, genome_to_cds_map = gene_to_cds_info[gene_name]
        df.at[ridx, CDS_SEQUENCE] = _to_str_seq(cds_seq)

        # External mRNA selection is computed here for future use (not stored as full sequence by default)
        mrna_built = _build_spliced_mrna_from_exons(pre_mrna, locus_info.exon_indices)
        _mrna_for_features = gene_to_mrna_real.get(gene_name, mrna_built)
        _ = _mrna_for_features  # keep explicit, prevents "unused" confusion

        # Local CDS context if within CDS and mappable
        if (
            locus_info.gene_start <= genome_corrected_index <= locus_info.gene_end
            and genome_corrected_index in genome_to_cds_map
        ):
            df.at[ridx, IN_CODING_REGION] = True
            cds_idx = genome_to_cds_map[genome_corrected_index]

            for fs in flank_sizes_cds:
                start = max(0, cds_idx - fs)
                end = min(len(cds_seq), cds_idx + len(sense) + fs)
                df.at[ridx, f"local_coding_region_around_ASO_{fs}"] = _to_str_seq(cds_seq[start:end])

    # Optional binary flags for modeling
    if add_region_is_local_flags:
        for fs in flank_sizes_cds:
            local_col = f"local_coding_region_around_ASO_{fs}"
            df[f"region_is_local_{fs}"] = df[local_col].apply(
                lambda x: isinstance(x, str) and x.strip() != ""
            ).astype(int)

    # Optional filter: keep only rows where site was found
    if filter_not_found:
        df = df[df[SENSE_START] != -1].copy()

    return df
