"""Inhibition vs off-target burden: a binned trend line per feature, with on-target as a control.

Every off-target feature bends efficacy down; the on-target control (binding to the INTENDED
target) stays flat -- so the effect is off-target-specific, not a general binding-strength thing.
y = within-cohort residual inhibition (oligo - its cohort mean). Marker size scales with the bin's
sample size, so sparse high-burden bins (where the trend gets noisy) read as low-confidence."""
import numpy as np
import matplotlib.pyplot as plt
import _feat as f

f.setup_style()
df = f.load()
rng = np.random.default_rng(0)

feats = list(f.BURDEN) + list(f.CONTROL)
fig, ax = plt.subplots(2, 4, figsize=(16, 7.6), sharey=True)
ax = ax.ravel()
for a, name in zip(ax, feats):
    control = name in f.CONTROL
    col = f.BLUE if control else f.ACCENT
    b = f.burden_log10(df, name)
    r = df['inh_resid'].to_numpy('float64')
    idx = rng.choice(len(b), min(5000, len(b)), replace=False)
    a.scatter(b[idx], r[idx], s=4, alpha=.06, color=f.GREY, edgecolors='none', rasterized=True)
    bl = f.binned_resid(df, name, nbins=10)
    a.fill_between(bl.x, bl['mean'] - bl.se, bl['mean'] + bl.se, color=col, alpha=.20, zorder=2)
    a.plot(bl.x, bl['mean'], '-', color=col, lw=1.8, zorder=3)
    a.scatter(bl.x, bl['mean'], s=18 + 130 * np.sqrt(bl.n / bl.n.max()), color=col,
              zorder=4, edgecolors='white', linewidths=.5)
    for _, row in bl.iterrows():
        if row.n < 0.04 * bl.n.max():
            a.annotate(f"n={int(row.n)}", (row.x, row['mean']), fontsize=6.5, color='#777',
                       ha='center', va='bottom', xytext=(0, 4), textcoords='offset points')
    a.axhline(0, color='#888', lw=.8); a.set_ylim(-22, 12)
    a.set_title(('● ' if control else '') + name, fontsize=10.5, color=(f.BLUE if control else '#222'))
    a.text(0, 1.015, ('CONTROL — ' if control else '') + f.BURDEN_DESC[name], transform=a.transAxes,
           fontsize=7.3, color='#777', va='bottom')
for a in ax[len(feats):]:
    a.set_visible(False)
for i in (0, 4):
    ax[i].set_ylabel('within-cohort residual\ninhibition (%)')
fig.text(.5, .01, 'off-target burden  (log₁₀ of the score; 0 = none)   ·   marker size ∝ bin sample size',
         ha='center', fontsize=9.5, color='#555')
fig.suptitle('Off-target binding lowers efficacy — on-target binding (blue control) does not',
             x=.012, ha='left', fontweight='bold', fontsize=14)
plt.tight_layout(rect=(0, .03, 1, .96))
f.save(fig, "inhibition_vs_burden")
