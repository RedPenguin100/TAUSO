"""Panel c - what a shortlist built by each method actually delivers.

A correlation says how well a method orders a whole screen; a user only ever tests the top of that
order. Per screen we take each method's top 5% and record the knockdown those ASOs really achieved,
against a random pick and against the best any selection could have done.

All violins share one density-to-width scale rather than matplotlib's per-violin normalisation, so
widths are comparable between methods: every method is scored on the same screens and each density
integrates to 1, so at any knockdown level the widths are true relative densities.
"""

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from _shared import ACCENT, BLUE, CID, GREY, INHIBITION, INK, MIN_ROWS, TEAL, save, scores, style
from scipy.stats import gaussian_kde

ORDER = ["random", "OligoWalk·intra", "OligoAI", "TAUSO", "best possible"]
COLOUR = {"random": GREY, "OligoWalk·intra": TEAL, "OligoAI": BLUE, "TAUSO": ACCENT, "best possible": INK}
SHORTLIST = 0.05


def shortlist_knockdown(frame, method):
    """Per screen, the mean measured knockdown of the top 5% this method would have picked."""
    out = []
    for _, screen in frame.groupby(CID):
        knockdowns = screen[INHIBITION].to_numpy(float)
        if len(screen) < MIN_ROWS or screen[INHIBITION].nunique() <= 1:
            continue
        k = max(1, int(round(SHORTLIST * len(screen))))
        if method == "random":
            out.append(knockdowns.mean())
        elif method == "best possible":
            out.append(np.sort(knockdowns)[-k:].mean())
        else:
            score = screen[method]
            if score.isna().any():
                continue
            out.append(knockdowns[np.argsort(-score.to_numpy(float))[:k]].mean())
    return np.array(out)


def main():
    frame = scores()
    values = {name: shortlist_knockdown(frame, name) for name in ORDER}
    values = {name: v for name, v in values.items() if len(v)}

    style()
    figure, ax = plt.subplots(figsize=(6.6, 4.8))
    names = list(values)
    grid = np.linspace(0, 105, 256)
    densities = {name: gaussian_kde(values[name])(grid) for name in names}
    tallest = max(d.max() for d in densities.values())

    for i, name in enumerate(names):
        inside = (grid >= values[name].min()) & (grid <= values[name].max())
        width = densities[name] / tallest * 0.42
        ax.fill_betweenx(grid[inside], i - width[inside], i + width[inside], facecolor=COLOUR[name],
                         alpha=0.5, edgecolor=COLOUR[name], linewidth=0.8, zorder=3)
        low, median, high = np.percentile(values[name], [25, 50, 75])
        ax.plot([i, i], [low, high], color="#222", lw=5, solid_capstyle="butt", zorder=4)
        ax.scatter([i], [median], color="white", s=14, edgecolor="#222", lw=0.5, zorder=5)
        if len(values[name]) < len(values["TAUSO"]):
            ax.text(i, 1.5, f"n={len(values[name])}", ha="center", va="bottom", fontsize=6.5, color="#555")

    halo = [pe.withStroke(linewidth=1.8, foreground="white")]
    here = names.index("TAUSO")
    low, median, high = np.percentile(values["TAUSO"], [25, 50, 75])
    ax.axhline(median, color=ACCENT, lw=1.0, alpha=0.40, zorder=2)
    ax.text(here + 0.13, median, f"{median:.1f}", color=INK, fontsize=8.5, va="center", ha="left",
            fontweight="bold", zorder=7, path_effects=halo)
    for quartile in (low, high):
        ax.axhline(quartile, color=ACCENT, ls="--", lw=0.7, alpha=0.28, zorder=2)
        ax.text(here + 0.13, quartile, f"{quartile:.1f}", color="#9aa3b0", fontsize=6.5, va="center",
                ha="left", zorder=7, path_effects=halo)

    ax.set_xticks(np.arange(len(names)))
    ax.set_xticklabels([n.replace(" ", "\n").replace("·", "\n·") for n in names], fontsize=8.5)
    ax.set_xlim(-0.6, len(names) - 0.4)
    ax.set_ylim(0, 105)
    ax.set_ylabel("top-5% knockdown (%), per screen", fontsize=10)
    ax.set_title("Top of the list: TAUSO's shortlist has a higher floor", fontsize=10.5,
                 fontweight="bold", loc="left", pad=6)
    save(figure, "panel_c_shortlist",
         caption="**c** Per screen, the measured knockdown of each method's top-5%-ranked ASOs. "
                 "Violin = distribution over screens on a shared density scale, bar = IQR, dot = median.")
    for name in names:
        print(f"  {name:18} median {np.median(values[name]):.1f}%  n={len(values[name])}")


if __name__ == "__main__":
    main()
