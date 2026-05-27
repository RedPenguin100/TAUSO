"""Cell-line name normalization for matching the TTDB half-life dataset."""

# Dataset spellings that differ from TTDB only by casing or punctuation. Names
# that already match TTDB (e.g. HepG2, HEK293T) are not listed and pass through.
# An entry is added only when the SAME line exists in TTDB under WT conditions
# (the same-line-or-nothing policy of CELL_LINE_TO_DEPMAP_PROXY_DICT); any cell
# line without a same-line TTDB match passes through unchanged and resolves to
# the gene-level estimate rather than a cross-lineage proxy.
_TTDB_CELL_NAME_FIXES = {
    "HeLa": "Hela",
    "MCF7": "MCF-7",
    "K-562": "K562",
}


def cell_name_to_ttdb(cell_name):
    """Return the TTDB spelling for a dataset cell-line name.

    Applies a known casing/punctuation fix when one exists; otherwise returns the
    name unchanged so the caller's lookup matches TTDB directly or falls back to
    the gene-level estimate.
    """
    return _TTDB_CELL_NAME_FIXES.get(cell_name, cell_name)
