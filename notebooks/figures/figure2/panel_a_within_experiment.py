"""Panel a - how well each method orders the ASOs inside one screen.

A screen (custom_id) is one gene, one cell line, one dose, assayed together, so ordering within it
is the question a design tool is actually asked. The bar is the median over screens and the whisker
a bootstrap interval over screens, not over ASOs: a method is being judged on how often it ranks a
screen well, not on how many ASOs it happens to cover.
"""

import matplotlib.pyplot as plt
from _shared import CID, median_ci, methods, per_screen_spearman, save, scores, spearman_bars, style


def main():
    frame = scores()
    medians, intervals = {}, {}
    for method in methods(frame):
        median, low, high = median_ci(per_screen_spearman(frame, method, CID))
        medians[method], intervals[method] = median, (low, high)

    style()
    figure, ax = plt.subplots(figsize=(6.4, 4.4))
    spearman_bars(ax, medians, intervals, "Within-experiment ranking (per screen)")
    save(figure, "panel_a_within_experiment",
         caption="**a** Median within-experiment Spearman between each method's score and measured "
                 "knockdown, over screens of at least 3 ASOs; whisker = 95% bootstrap CI over screens.")
    for method, value in sorted(medians.items(), key=lambda kv: -kv[1]):
        print(f"  {method:18} {value:.3f}  [{intervals[method][0]:.3f}, {intervals[method][1]:.3f}]")


if __name__ == "__main__":
    main()
