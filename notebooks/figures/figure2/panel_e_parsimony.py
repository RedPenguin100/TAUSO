"""Panel e - how much of the ranking survives when features are taken away.

Features are ranked by gain and a model retrained on the top K for decreasing K. The curve says how
concentrated the signal is: where it stays flat, the features being dropped were carrying nothing
the rest do not already carry.

The sweep is expensive, so it is precomputed into data/parsimony.csv rather than run here.
"""

import matplotlib.pyplot as plt
from _shared import ACCENT, BLUE, CID, median_ci, per_screen_spearman, save, scores, style, table
from matplotlib.ticker import FixedFormatter, FixedLocator, NullFormatter

TICKS = [20, 30, 50, 100, 200, 500]
ANNOTATE = (20, 30, 40, 50)


def main():
    curve = table("parsimony.csv").sort_values("K")
    baseline = median_ci(per_screen_spearman(scores(), "OligoAI", CID))[0]
    at = dict(zip(curve["K"], curve["exp_med"]))

    style()
    figure, ax = plt.subplots(figsize=(6.2, 4.4))
    ax.axhline(baseline, color=BLUE, lw=1.3, ls=(0, (5, 4)), zorder=2,
               label=f"strongest published baseline ({baseline:.2f})")
    ax.plot(curve["K"], curve["exp_med"], "-o", color=ACCENT, lw=1.6, ms=3.5, zorder=3, label="TAUSO")
    ax.scatter([curve["K"].max()], [at[curve["K"].max()]], marker="*", s=170, color=ACCENT,
               edgecolor="#222", lw=0.6, zorder=4)
    for k in ANNOTATE:
        if k in at:
            ax.annotate(f"{at[k]:.2f}", (k, at[k]), textcoords="offset points", xytext=(0, -12),
                        ha="center", fontsize=7.5, color=ACCENT, fontweight="bold")

    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator(TICKS))
    ax.xaxis.set_major_formatter(FixedFormatter([str(k) for k in TICKS]))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlim(17, curve["K"].max() * 1.15)
    ax.set_xlabel("number of features (top K by gain)", fontsize=10)
    ax.set_ylabel("within-screen median Spearman", fontsize=10)
    ax.set_title("Parsimony: few features suffice", fontsize=10.5, fontweight="bold", loc="left", pad=6)
    ax.legend(frameon=False, fontsize=8.5, loc="lower right")
    save(figure, "panel_e_parsimony",
         caption="**e** Within-screen median Spearman of models retrained on the top K gain-ranked "
                 "features; star = the shipped model.")
    print(f"  K range {int(curve['K'].min())}-{int(curve['K'].max())}, "
          f"plateau {curve['exp_med'].max():.3f}, baseline {baseline:.3f}")


if __name__ == "__main__":
    main()
