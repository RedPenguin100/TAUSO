"""Fig 2 — TAUSO vs. the field (held-out performance). Assembles panels a-f.

Data prep lives in data.py; each panel's drawing lives in panel_<x>.py. This file only lays out the
grid, calls each panel, and writes the figure. Heavy inputs (parsimony curve, screening-effort folds)
are precomputed by models/parsimony_gain.py and models/screening_effort.py.

    python figure.py
"""
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))                       # flat imports: data, panel_*, _helpers
sys.path.insert(0, str(HERE.parents[1]))            # notebooks/graphs, for consts (save/style)
from consts import save, style   # noqa: E402
import data as fig2data          # noqa: E402
import panel_a, panel_b, panel_c, panel_d, panel_e, panel_f   # noqa: E402


def _caption(D):
    eff = {m: D.effort_stat[m][0] for m in ("TAUSO", "OligoAI", "random")}
    return (
        "# Fig 2 — TAUSO vs. the field (held-out performance)\n\n"
        f"OligoAI's held-out test split (identical partition; Methods), strict 5-10-5 2′-MOE / 3-10-3 "
        f"cEt gapmers, rows OligoAI also scored (n = {D.n_aso:,} ASOs / {D.n_cid} patents / "
        f"{D.n_gxc} gene×cell cohorts). Both models are evaluated on identical rows unseen in training.\n\n"
        "**(a, b)** Median within-experiment and within-cohort Spearman vs measured knockdown; TAUSO "
        "leads both. **(c)** Per experiment, the actual knockdown of each method's top-5%-ranked ASOs "
        "(violin = distribution across experiments, bar = IQR, dot = median), from the random floor to "
        "the best possible. **(d)** Per-chemistry within-experiment Spearman; TAUSO's margin over "
        "OligoAI is largest on cEt. **(e)** Parsimony: retraining on the top-K gain-ranked features, "
        f"held-out Spearman holds at the shipped level down to ~80 features and beats OligoAI "
        f"({D.oai_em:.2f}) with ~30. **(f)** Screening-effort reduction vs random selection to match a "
        f"method's top-{D.effort_k} picks (geometric mean over screens, whisker = 95% bootstrap CI over "
        f"screens; fold is floored at 1×, so the random control sits just above it): TAUSO "
        f"{eff['TAUSO']:.1f}× vs OligoAI {eff['OligoAI']:.1f}× vs random {eff['random']:.1f}×. Draft."
    )


def main():
    D = fig2data.load()
    style()
    fig = plt.figure(figsize=(13.4, 14.6))
    gs = fig.add_gridspec(3, 2)
    panel_a.draw(fig.add_subplot(gs[0, 0]), D)
    panel_b.draw(fig.add_subplot(gs[0, 1]), D)
    panel_c.draw(fig.add_subplot(gs[1, 0]), D)
    panel_d.draw(fig.add_subplot(gs[1, 1]), D)
    panel_e.draw(fig.add_subplot(gs[2, 0]), D)
    panel_f.draw(fig.add_subplot(gs[2, 1]), D)
    fig.subplots_adjust(top=0.915, bottom=0.045, left=0.125, right=0.975, hspace=0.5, wspace=0.32)
    fig.suptitle("TAUSO vs. the field — held-out performance", fontsize=16, fontweight="bold",
                 x=0.012, y=0.975, ha="left")
    fig.text(0.012, 0.944, f"OligoAI split, strict 5-10-5 MOE / 3-10-3 cEt gapmers: "
             f"n = {D.n_aso:,} ASOs / {D.n_cid} patents / {D.n_gxc} gene×cell cohorts",
             fontsize=9.5, color="#666")
    save(fig, "article", "fig2_performance", caption=_caption(D))
    print(f"fig2 built | n={D.n_aso} TAUSO cid={D.cid['TAUSO']:.3f} gxc={D.gxc['TAUSO']:.3f} | "
          f"effort k={D.effort_k}")


if __name__ == "__main__":
    main()
