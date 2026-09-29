"""Panel f - how much wet-lab work each method saves.

To match the knockdown of a method's top five picks, how many ASOs would a random screen have to
test, keeping its own best five? The ratio is that method's saving, computed per screen and reported
as a geometric mean because the quantity is a ratio. The fold cannot fall below one, which is why
random sits a little above the 1x line rather than on it.
"""

import matplotlib.pyplot as plt
import numpy as np
from _shared import ACCENT, BLUE, GREY, save, style, table

ORDER = ["TAUSO", "OligoAI", "random"]
COLOUR = {"TAUSO": ACCENT, "OligoAI": BLUE, "random": GREY}


def geometric_mean_ci(folds, n_boot=4000):
    """Geometric mean of a method's per-screen folds, with a bootstrap interval over screens."""
    rng = np.random.default_rng(0)
    draws = np.exp(np.log(folds[rng.integers(0, len(folds), size=(n_boot, len(folds)))]).mean(axis=1))
    low, high = np.percentile(draws, [2.5, 97.5])
    return float(np.exp(np.log(folds).mean())), float(low), float(high)


def main():
    effort = table("screening_effort.csv")
    stats = {m: geometric_mean_ci(effort[m].to_numpy(float)) for m in ORDER}

    style()
    figure, ax = plt.subplots(figsize=(5.4, 4.4))
    x = np.arange(len(ORDER))
    means = [stats[m][0] for m in ORDER]
    ax.bar(x, means, 0.6, color=[COLOUR[m] for m in ORDER], edgecolor="white", zorder=3)
    error = np.array([[stats[m][0] - stats[m][1], stats[m][2] - stats[m][0]] for m in ORDER]).T
    ax.errorbar(x, means, yerr=error, fmt="none", ecolor="#444", elinewidth=1.1, capsize=4, zorder=4)
    ax.axhline(1.0, color="#888", lw=0.9, ls="--", zorder=2)

    tallest = max(stats[m][2] for m in ORDER)
    for i, method in enumerate(ORDER):
        ax.text(i, stats[method][2] + 0.03 * tallest, f"{stats[method][0]:.1f}x", ha="center",
                fontsize=11, fontweight="bold", color="#333")
    ax.set_xticks(x)
    ax.set_xticklabels(ORDER, fontsize=10)
    ax.set_ylim(0, tallest * 1.16)
    ax.set_ylabel("screening-effort reduction (x)\nto match the top-5 picks", fontsize=10)
    ax.set_title("Screening effort saved", fontsize=10.5, fontweight="bold", loc="left", pad=6)
    save(figure, "panel_f_screening_effort",
         caption="**f** Factor by which each method shrinks a random screen matching its top-5 picks "
                 f"(geometric mean over {len(effort)} screens; whisker = 95% bootstrap CI over screens).")
    for method in ORDER:
        mean, low, high = stats[method]
        print(f"  {method:8} {mean:.2f}x  [{low:.2f}, {high:.2f}]")


if __name__ == "__main__":
    main()
