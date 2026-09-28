"""Panel b - how well each method orders ASOs across a whole gene x cell-line cohort.

A cohort pools every screen run on one gene in one cell line, so it spans different doses, delivery
routes and gapmer architectures. Ranking across those is a harder question than ranking inside a
single screen, and a method that has only learned the idiosyncrasies of individual screens loses
its advantage here.
"""

import matplotlib.pyplot as plt
from _shared import COHORT, median_ci, methods, per_screen_spearman, save, scores, spearman_bars, style


def main():
    frame = scores()
    medians, intervals = {}, {}
    for method in methods(frame):
        median, low, high = median_ci(per_screen_spearman(frame, method, COHORT))
        medians[method], intervals[method] = median, (low, high)

    style()
    figure, ax = plt.subplots(figsize=(6.4, 4.4))
    spearman_bars(ax, medians, intervals, "Within-cohort ranking (per gene x cell line)")
    save(figure, "panel_b_within_cohort",
         caption="**b** Median within-cohort Spearman between each method's score and measured "
                 "knockdown, over cohorts of at least 3 ASOs; whisker = 95% bootstrap CI over cohorts.")
    for method, value in sorted(medians.items(), key=lambda kv: -kv[1]):
        print(f"  {method:18} {value:.3f}  [{intervals[method][0]:.3f}, {intervals[method][1]:.3f}]")


if __name__ == "__main__":
    main()
