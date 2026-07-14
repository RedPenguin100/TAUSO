"""Recompute fold features for splice-junction ASOs against the spliced canonical CDS.

A splice-junction ASO targets an exon-exon boundary, so its target is contiguous only in the
spliced transcript and is split by an intron in the gene-span pre-mRNA (``full_mrna``). It
therefore fails to locate on the pre-mRNA (``structure_sense_start == -1``) and its fold features
(``fold_mfe_*`` / ``fold_access_*``) come out NaN. Locating the ASO on the canonical CDS instead
yields a CDS-frame start, and the fold populate functions fold it against that spliced template.
"""

import numpy as np

from ..data.consts import CANONICAL_GENE_NAME, STRUCTURE_SENSE_LENGTH, STRUCTURE_SENSE_START
from ..genome.LocusInfo import CanonicalCds
from ..util import get_antisense_rna
from .populate_fold import (
    DEFAULT_SENSE_CONFIGURATION,
    SEQUENCE_SOURCE_CANONICAL_CDS,
    populate_mfe_features,
    populate_sense_accessibility_multi,
)
from .populate_structure import _match_positions

FOLD_SOURCE_COL = "fold_source"
_T_TO_U = str.maketrans("tT", "uU")


def _cds_target_bytes(cds_sequence):
    """CDS as the same uppercase RNA-alphabet byte-string the pre-mRNA locator uses as its target."""
    return cds_sequence.upper().translate(_T_TO_U).encode()


def locate_on_cds(sequences, genes, gene_to_data):
    """CDS-frame start of each ASO's target in its gene's canonical CDS, or -1 if not found there.

    Mirrors the pre-mRNA locator (search the antisense-RNA of the ASO with ``_match_positions``),
    but searches ``CanonicalCds.from_locus(locus).sequence`` instead of ``full_mrna``. The CDS is
    built once per gene. Genes absent from ``gene_to_data`` or with an empty canonical CDS yield -1.
    """
    sequences = list(sequences)
    genes = list(genes)
    starts = np.full(len(sequences), -1, dtype=np.int64)

    rows_by_gene = {}
    for i, gene in enumerate(genes):
        rows_by_gene.setdefault(gene, []).append(i)

    for gene, rows in rows_by_gene.items():
        locus = gene_to_data.get(gene)
        if locus is None:
            continue
        cds_sequence = CanonicalCds.from_locus(locus).sequence
        if not cds_sequence:
            continue
        target = _cds_target_bytes(cds_sequence)
        antisenses = [get_antisense_rna(sequences[i]).encode() for i in rows]
        positions = _match_positions(target, antisenses)
        for i, pos in zip(rows, positions):
            starts[i] = int(pos)
    return starts


def recompute_fold_on_cds(
    df,
    gene_to_data,
    sequence_col,
    n_jobs=1,
    mfe_settings=None,
    access_configs=None,
    with_mfe=True,
    with_access=True,
):
    """Fold the ASOs in ``df`` against their canonical CDS.

    Relocates each ASO on the CDS (``locate_on_cds``), writes the CDS-frame
    ``STRUCTURE_SENSE_START`` / ``STRUCTURE_SENSE_LENGTH`` and a ``FOLD_SOURCE_COL`` of
    ``canonical_cds``, then runs the fold populate functions with that per-row source. Rows that do
    not locate on the CDS keep ``structure_sense_start == -1`` and fall through to NaN fold features
    (the populate functions skip them), so a genuinely non-CDS target is left blank rather than wrong.

    Returns ``(df, fold_feature_names)``; ``df`` carries the fold columns plus the CDS start it used.
    """
    df = df.copy()
    df[STRUCTURE_SENSE_START] = locate_on_cds(
        df[sequence_col].to_numpy(), df[CANONICAL_GENE_NAME].to_numpy(), gene_to_data
    )
    df[STRUCTURE_SENSE_LENGTH] = df[sequence_col].str.len().to_numpy()
    df[FOLD_SOURCE_COL] = SEQUENCE_SOURCE_CANONICAL_CDS

    fold_feature_names = []
    if with_mfe:
        df, mfe_names = populate_mfe_features(
            df, gene_to_data, n_jobs=n_jobs, settings=mfe_settings, sequence_source_col=FOLD_SOURCE_COL
        )
        fold_feature_names += mfe_names
    if with_access:
        df, access_names = populate_sense_accessibility_multi(
            df,
            gene_to_data,
            access_configs if access_configs is not None else DEFAULT_SENSE_CONFIGURATION,
            n_jobs=n_jobs,
            sequence_source_col=FOLD_SOURCE_COL,
        )
        fold_feature_names += access_names
    return df, fold_feature_names
