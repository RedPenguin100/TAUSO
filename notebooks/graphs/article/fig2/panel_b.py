"""Fig 2b — within-cohort ranking (per gene x cell-line)."""
from _helpers import panel_label, spearman_bar


def draw(ax, D):
    spearman_bar(ax, D.gxc, "Within-cohort ranking (per gene × cell-line)")
    panel_label(ax, "b")
