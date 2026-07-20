"""Fig 2e — parsimony: held-out within-experiment Spearman vs number of top-K features."""
from _helpers import ACCENT, BLUE, INK, TEAL, droplabels, logx_features, panel_label

LOWK = [20, 30, 40, 50]


def draw(ax, D):
    pt = D.parsimony
    ptmap = dict(zip(pt["K"], pt["exp_med"]))
    kmax = int(pt["K"].max())
    ax.axhline(D.oai_em, color=BLUE, lw=1.3, ls=(0, (5, 4)), zorder=2, label=f"OligoAI ({D.oai_em:.2f})")
    ax.axhline(D.ow_em, color=TEAL, lw=1.3, ls=(0, (5, 4)), zorder=2, label=f"{D.ow_lbl} ({D.ow_em:.2f})")
    ax.plot(pt["K"], pt["exp_med"], color=ACCENT, lw=2.0, marker="o", ms=2.8, zorder=4, label="TAUSO")
    ax.scatter([kmax], [ptmap[kmax]], s=150, marker="*", color=ACCENT, edgecolor=INK, lw=0.8, zorder=5)
    droplabels(ax, LOWK, lambda k: ptmap[k], ACCENT)
    logx_features(ax, kmax)
    ax.set_ylabel("within-experiment Spearman", fontsize=10)
    ax.set_title("Parsimony — few features suffice", fontsize=10.5, fontweight="bold", loc="left", pad=6)
    ax.legend(frameon=False, fontsize=8.5, loc="lower right")
    panel_label(ax, "e")
