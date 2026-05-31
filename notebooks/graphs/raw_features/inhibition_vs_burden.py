"""Inhibition vs off-target burden: faint per-oligo scatter + a binned trend line.

Shows the continuous relationship (binned, not percentile-split) and what burden magnitudes
correspond to what efficacy drop. y is within-cohort residual inhibition (oligo minus its cohort
mean), so the experiment is controlled."""
import numpy as np
import matplotlib.pyplot as plt
import _feat as f

f.setup_style()
df = f.load()
rng = np.random.default_rng(0)

feats = ['rRNA (total)', 'off-target general']
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for a, name in zip(ax, feats):
    b = f.burden_log10(df, name)
    r = df['inh_resid'].to_numpy('float64')
    idx = rng.choice(len(b), min(7000, len(b)), replace=False)
    a.scatter(b[idx], r[idx], s=5, alpha=.07, color=f.GREY, edgecolors='none', rasterized=True)
    bl = f.binned_resid(df, name, nbins=10)
    a.fill_between(bl.x, bl['mean'] - bl.se, bl['mean'] + bl.se, color=f.ACCENT, alpha=.22, zorder=2)
    a.plot(bl.x, bl['mean'], '-o', color=f.ACCENT, lw=2.4, ms=6, zorder=3)
    a.axhline(0, color='#888', lw=.8)
    a.set_ylim(-22, 12)
    a.set_title(name, pad=24)
    a.text(0, 1.012, f.BURDEN_DESC[name], transform=a.transAxes, fontsize=8.5, color='#666', va='bottom')
    a.set_xlabel('off-target burden  (log₁₀ of the score; 0 = none)')
ax[0].set_ylabel('within-cohort residual inhibition (%)\n(oligo − its cohort mean)')
fig.suptitle('More off-target binding → lower efficacy', x=.012, ha='left',
             fontweight='bold', fontsize=14, y=1.04)
plt.tight_layout()
f.save(fig, "inhibition_vs_burden")
