"""The expression tables the features read: how they are built, and where they are kept.

`general` holds the per-gene mean across the DepMap cohort, which the general off-target
features rank genes by. `cohort` holds the per-cell-line tables, gene level and transcript
level. Building them belongs here rather than beside the features that read them, so that
the calculator's cache only has to remember what it has already loaded.
"""
