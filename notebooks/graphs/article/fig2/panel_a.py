"""Fig 2a — within-experiment ranking (per patent / custom_id)."""
from _helpers import panel_label, spearman_bar


def draw(ax, D):
    spearman_bar(ax, D.cid, "Within-experiment ranking (per patent / custom_id)")
    panel_label(ax, "a")
