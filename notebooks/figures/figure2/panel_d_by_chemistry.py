"""Panel d - the same within-screen ranking, split by gapmer architecture.

The two canonical architectures differ in wing chemistry and in length, and a method can be good at
one and poor at the other. Splitting the metric says where an advantage actually comes from instead
of averaging the two into a single number.
"""

import matplotlib.pyplot as plt
import numpy as np
from _shared import ACCENT, BLUE, CID, median_ci, per_screen_spearman, save, scores, style

CHEMISTRIES = {"2'-MOE": "MMMMMddddddddddMMMMM", "cEt": "CCCddddddddddCCC"}
COMPARED = (("TAUSO", ACCENT), ("OligoAI", BLUE))


def main():
    frame = scores()
    stats = {}
    for name, pattern in CHEMISTRIES.items():
        part = frame[frame["chemical_pattern"] == pattern]
        stats[name] = {m: median_ci(per_screen_spearman(part, m, CID)) for m, _ in COMPARED}

    style()
    figure, ax = plt.subplots(figsize=(5.8, 4.4))
    x = np.arange(len(CHEMISTRIES))
    width = 0.36
    tallest = 0.0
    for offset, (method, colour) in enumerate(COMPARED):
        medians = np.array([stats[c][method][0] for c in CHEMISTRIES])
        bounds = [stats[c][method][1:] for c in CHEMISTRIES]
        error = np.array([[m - lo, hi - m] for m, (lo, hi) in zip(medians, bounds)]).T
        ax.bar(x + (offset - 0.5) * width, medians, width, color=colour, edgecolor="white", zorder=3,
               label=method, yerr=error, capsize=3,
               error_kw=dict(ecolor="#444", elinewidth=1.0, zorder=4))
        for i, (median, (_low, high)) in enumerate(zip(medians, bounds)):
            ax.text(x[i] + (offset - 0.5) * width, high + 0.012, f"{median:.2f}", ha="center",
                    fontsize=8.5, color="#333")
            tallest = max(tallest, high)

    ax.set_xticks(x)
    ax.set_xticklabels(list(CHEMISTRIES), fontsize=11)
    ax.set_ylim(0, tallest * 1.14)
    ax.set_ylabel("within-screen median Spearman", fontsize=10)
    ax.set_title("Per-chemistry ranking", fontsize=11, fontweight="bold", loc="left", pad=6)
    ax.legend(frameon=False, fontsize=9.5, loc="upper right")
    save(figure, "panel_d_by_chemistry",
         caption="**d** Within-screen median Spearman split by gapmer architecture; "
                 "whisker = 95% bootstrap CI over screens.")
    for name in CHEMISTRIES:
        gap = stats[name]["TAUSO"][0] - stats[name]["OligoAI"][0]
        print(f"  {name:8} TAUSO {stats[name]['TAUSO'][0]:.3f}  OligoAI {stats[name]['OligoAI'][0]:.3f}"
              f"  margin {gap:+.3f}")


if __name__ == "__main__":
    main()
