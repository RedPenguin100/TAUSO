"""Panel-c alternative: survival curve of the top-5% shortlist quality. For each method, per experiment take
the mean knockdown of its top-5%-ranked ASOs; plot the fraction of experiments whose shortlist averages >= x,
as x sweeps 0..100. Higher/flatter-until-late = a higher floor (more consistent). Preview for fig 2 panel c.
Run: python fig2c_survival.py
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "models"))
from consts import ACCENT, BLUE, GREY, INK
import _modeldata as M

INH, CID = M.INH, M.CID
S, present = M.load()


def shortlist(score):
    out = []
    for _, g in S.groupby(CID):
        n = len(g)
        if n >= 10 and g[INH].nunique() > 1:
            k = max(1, int(round(0.05 * n)))
            v = (g[INH].mean() if score == "_random" else
                 g.nlargest(k, INH)[INH].mean() if score == "_oracle" else
                 g.nlargest(k, score)[INH].mean())
            out.append(v)
    return np.array(out)


methods = [("best possible", "_oracle", INK, "--"), ("TAUSO", "tauso", ACCENT, "-"),
           ("OligoAI", "oligo_ai_score", BLUE, "-"), ("Random", "_random", GREY, "-")]
xs = np.linspace(0, 100, 201)
fig, ax = plt.subplots(figsize=(6.2, 5.0))
vals = {}
for name, score, c, ls in methods:
    a = shortlist(score); vals[name] = a
    ax.plot(xs, [100 * np.mean(a >= x) for x in xs], color=c, lw=2.4, ls=ls, label=name, zorder=3)

# annotate the floor gap at 50%
for name, c in [("TAUSO", ACCENT), ("OligoAI", BLUE)]:
    f50 = 100 * np.mean(vals[name] >= 50)
    ax.plot([50], [f50], "o", color=c, ms=6, zorder=4)
ax.axvline(50, color="#bbb", lw=0.8, ls=":", zorder=1)
gap = 100 * (np.mean(vals["TAUSO"] >= 50) - np.mean(vals["OligoAI"] >= 50))
ax.annotate(f"at ≥50%: TAUSO {100*np.mean(vals['TAUSO']>=50):.0f}% vs OligoAI {100*np.mean(vals['OligoAI']>=50):.0f}%\n(TAUSO duds half as often)",
            xy=(50, 100 * np.mean(vals["OligoAI"] >= 50)), xytext=(8, 22), fontsize=8.5, color="#333")

ax.set_xlabel("shortlist mean knockdown ≥ x  (%)")
ax.set_ylabel("% of experiments")
ax.set_title("Top-5% shortlist: TAUSO holds a higher floor\n(curves converge at the top = tied ceiling; TAUSO wins the middle)", fontsize=11)
ax.set_xlim(0, 100); ax.set_ylim(0, 101)
ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
ax.grid(alpha=0.25)
fig.tight_layout()
out = Path(__file__).resolve().parents[1] / "out/article/fig2c_survival.png"
fig.savefig(out, dpi=145, bbox_inches="tight")
print("wrote", out)
for name in ["TAUSO", "OligoAI"]:
    a = vals[name]
    print(f"  {name}: ≥80 {100*np.mean(a>=80):.0f}% | ≥60 {100*np.mean(a>=60):.0f}% | ≥50 {100*np.mean(a>=50):.0f}% | ≥40 {100*np.mean(a>=40):.0f}%")
