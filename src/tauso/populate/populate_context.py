import logging
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd

from ..data.consts import CANONICAL_GENE_NAME, CELL_LINE_DEPMAP, TRANSFECTION_RAW

logger = logging.getLogger(__name__)


# Proteins a PS-ASO has to get past in the cell it is dosed into, taken per gene. EGFR,
# stabilin-1 and stabilin-2 are endocytic receptors for phosphorothioate oligonucleotides;
# TREX1 is the 3'-exonuclease that degrades them. A gene total is the measure here because
# these features stand for how much of the protein a cell line carries, which does not turn on
# which isoform the message sits in.
_SPECIAL_GENES: Dict[str, str] = {
    "expr_egfr": "EGFR",
    "expr_stab1": "STAB1",
    "expr_stab2": "STAB2",
    "expr_trex1": "TREX1",
}

TARGET_EXPRESSION_FEATURE_NAMES: List[str] = ["expr_target"]

# How concentrated the target's expression is: the share of its TPM carried by its single most
# abundant transcript, from 0 to 1. A gene at 100 TPM in one isoform and one at 100 TPM split
# across five present the ASO with different amounts of the sequence it was designed against,
# and expr_target alone cannot tell them apart. A share rather than an amount, so it says
# something expr_target does not already carry. The dominant transcript is resolved per cell
# line rather than named: the target differs from row to row.
TARGET_TRANSCRIPT_FEATURE_NAMES: List[str] = ["expr_target_dom_fraction"]

# RNase H1 cleaves the ASO:RNA heteroduplex, and unlike the genes above it is taken per
# transcript: RNASEH1-201 carries 90% of the gene, so the gene total and the canonical
# transcript say the same thing and the transcript says it about one species.
#
# A feature names every transcript that encodes the species it is about, and their expression is
# summed. Transcripts are named outright rather than resolved per cell line, so a column holds
# the same species everywhere. Ensembl transcript names are used rather than accessions because
# they say which gene and which isoform without a lookup. A name that falls out of a future
# annotation leaves the feature NaN and says so, rather than quietly switching isoform.
_SPECIAL_TRANSCRIPTS: Dict[str, Tuple[str, ...]] = {
    "expr_rnase_transcript": ("RNASEH1-201",),
}
SPECIAL_TRANSCRIPT_EXPRESSION_FEATURE_NAMES: List[str] = list(_SPECIAL_TRANSCRIPTS.keys())
SPECIAL_GENE_EXPRESSION_FEATURE_NAMES: List[str] = list(_SPECIAL_GENES.keys())
EXPRESSION_FEATURE_NAMES: List[str] = (
    TARGET_EXPRESSION_FEATURE_NAMES
    + TARGET_TRANSCRIPT_FEATURE_NAMES
    + SPECIAL_GENE_EXPRESSION_FEATURE_NAMES
    + SPECIAL_TRANSCRIPT_EXPRESSION_FEATURE_NAMES
)


def _build_expression_master(expression_dict: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Convert from expression_dict (depmap_id : DataFrame(gene / expr)
    To expression_master (DataFrame (gene / expr / depmap_id))
    """
    dfs = []
    for depmap_id, t_df in expression_dict.items():
        temp = t_df[["Gene", "expression_norm"]].copy()
        temp[CELL_LINE_DEPMAP] = depmap_id
        dfs.append(temp)
    master = pd.concat(dfs, ignore_index=True)
    return master.drop_duplicates(subset=[CELL_LINE_DEPMAP, "Gene"])


def _drop_existing(df: pd.DataFrame, cols: Sequence[str]) -> pd.DataFrame:
    existing_cols = [c for c in cols if c in df.columns]
    if existing_cols:
        logger.warning("Dropping existing columns: %s", existing_cols)
        return df.drop(columns=existing_cols)
    return df


def populate_target_expression(
    df: pd.DataFrame, expression_dict: Dict[str, pd.DataFrame]
) -> Tuple[pd.DataFrame, List[str]]:
    """Merges the ASO's target gene expression (per cell line × gene) into df."""
    _drop_existing(df, TARGET_EXPRESSION_FEATURE_NAMES)

    expression_master = _build_expression_master(expression_dict)

    enhanced_df = df.merge(
        expression_master,
        left_on=[CELL_LINE_DEPMAP, CANONICAL_GENE_NAME],
        right_on=[CELL_LINE_DEPMAP, "Gene"],
        how="left",
    )
    enhanced_df = enhanced_df.rename(columns={"expression_norm": "expr_target"})
    enhanced_df = enhanced_df.drop(columns=["Gene"])

    return enhanced_df, TARGET_EXPRESSION_FEATURE_NAMES


def populate_target_dominant_transcript(
    df: pd.DataFrame, transcript_dict: Dict[str, pd.DataFrame]
) -> Tuple[pd.DataFrame, List[str]]:
    """Merges the target gene's dominant-transcript share (per cell line x gene) into df.

    The share is the most abundant transcript's TPM over the gene's total, so 1.0 is a gene
    expressed as a single species and 0.2 one spread across several.
    """
    _drop_existing(df, TARGET_TRANSCRIPT_FEATURE_NAMES)

    frames = []
    for depmap_id, t_df in transcript_dict.items():
        grouped = t_df.groupby("Gene")["expression_TPM"]
        best = grouped.max().rename("top_tpm").to_frame()
        best["gene_tpm"] = grouped.sum()
        best = best.reset_index().assign(**{CELL_LINE_DEPMAP: depmap_id})
        frames.append(best)

    if not frames:
        logger.warning("No transcript expression found for the target genes")
        df[TARGET_TRANSCRIPT_FEATURE_NAMES[0]] = np.nan
        return df, TARGET_TRANSCRIPT_FEATURE_NAMES

    master = pd.concat(frames, ignore_index=True)
    # A gene with no expression at all has no dominant share to speak of.
    master["expr_target_dom_fraction"] = np.where(
        master["gene_tpm"] > 0, master["top_tpm"] / master["gene_tpm"], np.nan
    )
    master = master.drop(columns=["top_tpm", "gene_tpm"])

    merged = df.merge(
        master,
        left_on=[CELL_LINE_DEPMAP, CANONICAL_GENE_NAME],
        right_on=[CELL_LINE_DEPMAP, "Gene"],
        how="left",
    ).drop(columns=["Gene"])

    return merged, TARGET_TRANSCRIPT_FEATURE_NAMES


def _merge_fixed_gene(
    df: pd.DataFrame, expression_master: pd.DataFrame, gene_symbol: str, feat_name: str
) -> pd.DataFrame:
    gene_vals = expression_master[expression_master["Gene"] == gene_symbol]
    merged = df.merge(
        gene_vals[[CELL_LINE_DEPMAP, "expression_norm"]],
        on=CELL_LINE_DEPMAP,
        how="left",
        suffixes=("", f"_{feat_name}"),
    )
    return merged.rename(columns={"expression_norm": feat_name})


def populate_special_gene_expression(
    df: pd.DataFrame, expression_dict: Dict[str, pd.DataFrame]
) -> Tuple[pd.DataFrame, List[str]]:
    """Merges per-cell-line expression for the genes named in _SPECIAL_GENES."""
    _drop_existing(df, SPECIAL_GENE_EXPRESSION_FEATURE_NAMES)

    expression_master = _build_expression_master(expression_dict)

    for feat_name, gene_symbol in _SPECIAL_GENES.items():
        df = _merge_fixed_gene(df, expression_master, gene_symbol, feat_name)

    return df, SPECIAL_GENE_EXPRESSION_FEATURE_NAMES


def populate_special_transcript_expression(
    df: pd.DataFrame, expression_dict: Dict[str, pd.DataFrame]
) -> Tuple[pd.DataFrame, List[str]]:
    """Merges per-cell-line expression of each named special transcript.

    The transcript-level twin of populate_special_gene_expression. Where that one takes a
    gene's expression, this takes one named isoform of it: a gene whose expression sits in a
    single isoform and one split across four share a gene-level value but hold different
    amounts of any one species. Ids are matched unversioned.
    """
    _drop_existing(df, SPECIAL_TRANSCRIPT_EXPRESSION_FEATURE_NAMES)

    rows = []
    for depmap_id, t_df in expression_dict.items():
        indexed = t_df.set_index("TranscriptName")
        for feat_name, transcript_names in _SPECIAL_TRANSCRIPTS.items():
            present = [n for n in transcript_names if n in indexed.index]
            if not present:
                continue
            if len(present) < len(transcript_names):
                logger.warning(
                    "%s: %s absent for %s", feat_name, sorted(set(transcript_names) - set(present)), depmap_id
                )
            # Expression is log2(TPM + 1), so it is summed in TPM and converted back.
            total_tpm = float(indexed.loc[present, "expression_TPM"].sum())
            rows.append({"Gene": feat_name, "expression_norm": np.log2(total_tpm + 1.0), CELL_LINE_DEPMAP: depmap_id})

    if not rows:
        logger.warning("No expression found for the special transcripts")
        for feat_name in SPECIAL_TRANSCRIPT_EXPRESSION_FEATURE_NAMES:
            df[feat_name] = np.nan
        return df, SPECIAL_TRANSCRIPT_EXPRESSION_FEATURE_NAMES

    master = pd.DataFrame(rows)
    for feat_name in _SPECIAL_TRANSCRIPTS:
        df = _merge_fixed_gene(df, master, feat_name, feat_name)

    return df, SPECIAL_TRANSCRIPT_EXPRESSION_FEATURE_NAMES


_TRANSFECTION_LABEL_TO_COLUMN = {
    "Electroporation": "transfection_electroporation",
    "Gymnosis": "transfection_gymnosis",
    "Lipofection": "transfection_lipofection",
}


def populate_transfection(data):
    """
    One-hot encodes the raw transfection_raw label into three columns:
    transfection_electroporation / transfection_gymnosis / transfection_lipofection.
    Rows whose transfection_raw isn't one of those three (e.g. "Other", missing,
    or any unrecognized label) get NaN across all three -- the model treats them
    as missing rather than as confidently "not Electroporation, not Gymnosis,
    not Lipofection".

    Returns float64 columns because NaN cannot live in an int column.
    """
    features = list(_TRANSFECTION_LABEL_TO_COLUMN.values())
    labels = list(_TRANSFECTION_LABEL_TO_COLUMN.keys())

    method = data[TRANSFECTION_RAW]
    binary = pd.DataFrame(
        {col: (method == label).astype("float64") for label, col in _TRANSFECTION_LABEL_TO_COLUMN.items()},
        index=data.index,
    )
    unknown = ~method.isin(labels)
    binary.loc[unknown, features] = np.nan

    data = pd.concat([data, binary], axis=1)
    return data, features
