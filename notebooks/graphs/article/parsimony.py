"""Parsimony -- held-out within-experiment Spearman vs number of top-K features (gain-ranked).

Reads the sweep written by models/parsimony_gain.py and draws the curve against the competitor
reference lines, on the same held-out comparison subset Fig 2 uses.

    python parsimony.py
"""
import sys
from pathlib import Path

import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from consts import save, style, ACCENT, BLUE, GREY, INK   # noqa: E402

CAT = "article"
MODELS = Path(__file__).resolve().parents[1] / "models"
LOWK = [20, 30, 40, 50]
XTICKS = [20, 30, 50, 100, 200, 500]
TEAL = "#2A9D8F"


def main():
    d = pd.read_csv(MODELS / "parsimony_gain_test.csv").sort_values("K")
    refs = pd.read_csv(MODELS / "parsimony_gain_refs.csv").set_index("label")["exp_med"]
    K, em = d["K"].to_numpy(), d["exp_med"].to_numpy()
    at = dict(zip(d["K"], d["exp_med"]))
    kmax = int(d["K"].max())

    oai = float(refs.get("OligoAI", float("nan")))
    rule_lbl = refs.drop(index=["OligoAI"], errors="ignore").idxmax()   # best rule-based tool
    rule = float(refs[rule_lbl])

    style()
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.axhline(oai, color=BLUE, lw=1.3, ls=(0, (5, 4)), zorder=2, label=f"OligoAI ({oai:.2f})")
    ax.axhline(rule, color=TEAL, lw=1.3, ls=(0, (5, 4)), zorder=2, label=f"{rule_lbl} ({rule:.2f})")
    ax.plot(K, em, color=ACCENT, lw=2.0, marker="o", ms=3.2, zorder=4,
            label=f"TAUSO (K={kmax}: {at[kmax]:.2f})")
    ax.scatter([kmax], [at[kmax]], s=170, marker="*", color=ACCENT, edgecolor=INK, lw=0.8, zorder=5)
    ax.annotate(f"shipped model\n(all {kmax} features)", (kmax, at[kmax]), textcoords="offset points",
                xytext=(-10, -30), ha="right", fontsize=8, color=ACCENT, fontweight="bold")
    for k in LOWK:
        if k in at:
            ax.annotate(f"{at[k]:.2f}", (k, at[k]), textcoords="offset points", xytext=(0, -13),
                        ha="center", fontsize=8, color=ACCENT, fontweight="bold")

    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator(XTICKS))
    ax.xaxis.set_major_formatter(FixedFormatter([str(k) for k in XTICKS]))
    ax.set_xlim(17, kmax * 1.18)
    ax.set_xlabel("number of features (top-K, ranked by gain)", fontsize=10)
    ax.set_ylabel("within-experiment Spearman", fontsize=10)
    ax.set_title("Parsimony — held-out test", fontsize=12, fontweight="bold", loc="left", pad=8)
    ax.legend(frameon=False, fontsize=9, loc="lower right")

    cap = (
        "# Parsimony — held-out test\n\n"
        "Within-experiment (custom_id) median Spearman on the held-out test, restricted to the Fig 2 "
        "comparison subset (strict 5-10-5 2'-MOE / 3-10-3 cEt gapmers, rows OligoAI also scored). "
        "For each K the champion configuration is retrained on train+val using only the top-K features "
        "ranked by **gain** importance, then scored once on the test split. Dashed lines are the "
        f"competitors on the same rows: OligoAI ({oai:.2f}) and the strongest rule-based tool, "
        f"{rule_lbl} ({rule:.2f}). ★ marks the shipped model, which uses all {kmax} features "
        f"({at[kmax]:.2f})."
    )
    save(fig, CAT, "parsimony", caption=cap)
    print(f"parsimony built | K {K.min()}..{K.max()} ({len(K)} points) | "
          f"K={kmax}: {at[kmax]:.3f} | OligoAI {oai:.3f} | {rule_lbl} {rule:.3f}")


if __name__ == "__main__":
    main()
