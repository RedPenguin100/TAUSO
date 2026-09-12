"""Hybridization off-target / on-target site-feature configuration.

Single source of truth for the (top_n, cutoff) grid the calculator computes these
features on.
"""

from pyrisearch_tauso import Matrix

# RIsearch score cutoffs (a hit counts toward cutoff c when its score > c). Shared by
# the on-target, single-gene, and general/specific off-target hybridization features.
RISEARCH_SCORE_CUTOFFS = (800, 1000, 1200)

# Expression-ranked gene-universe sizes (head(top_n) of the transcriptome) for the
# general/specific off-target scores.
OFF_TARGET_TOP_NS = (50, 100, 200)

# Calibrated ASO-target search parameters and RT at body temperature, kcal/mol.
ASO_TARGET_MATRIX = Matrix.SU95_NO_GU
EXTENSION_PENALTY = 30  # dacal/mol per nucleotide
RT_KCAL_MOL = 0.616
