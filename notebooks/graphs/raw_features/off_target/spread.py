"""Spread/significance: per-cohort delta (high - regular) distribution, mean +/- 95% CI, Wilcoxon p."""
import numpy as np
import matplotlib.pyplot as plt
import _feat as f

f.setup_style()
df = f.load()

feats = ['rRNA (total)', 'rRNA 18S', 'RNase-H1 off-tgt', 'off-target general', 'off-target specific']
res = {name: f.cohort_delta(df, name) for name in feats}
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6), gridspec_kw={'width_ratios': [1.25, 1]})

c = res['rRNA (total)']; delta = c['delta']
ax[0].hist(delta, bins=40, color=f.GREY, alpha=.85, edgecolor='white', linewidth=.4)
ax[0].axvline(0, color='#333', lw=1.2)
ax[0].axvspan(c['mean'] - c['ci'], c['mean'] + c['ci'], color=f.ACCENT, alpha=.18)
ax[0].axvline(c['mean'], color=f.ACCENT, lw=2)
ax[0].annotate(f"mean {c['mean']:+.1f}\n95% CI [{c['mean']-c['ci']:+.1f}, {c['mean']+c['ci']:+.1f}]",
               (c['mean'], ax[0].get_ylim()[1] * .92), ha='right', va='top', fontsize=9, color=f.ACCENT)
ax[0].annotate(f"{100*(delta<0).mean():.0f}% of cohorts < 0", (.97, .80), xycoords='axes fraction',
               ha='right', fontsize=9, color='#555')
ax[0].set_xlabel('per-cohort delta: high − regular  (% inhibition)')
ax[0].set_ylabel('cohorts'); ax[0].set_title('rRNA (total): per-cohort deltas')

y = np.arange(len(feats))[::-1]
ax[1].axvline(0, color='#333', lw=1.2)
for yi, name in zip(y, feats):
    r = res[name]
    ax[1].errorbar(r['mean'], yi, xerr=r['ci'], fmt='o', ms=8, color=f.ACCENT,
                   ecolor=f.ACCENT, elinewidth=1.6, capsize=4)
    ax[1].annotate(f"p={r['p']:.0e}", (r['mean'], yi + 0.18), ha='center', fontsize=8, color='#555')
ax[1].set_yticks(y); ax[1].set_yticklabels(feats); ax[1].set_ylim(-0.6, len(feats) - 0.3)
ax[1].set_xlabel('mean delta ± 95% CI  (% inhibition)')
ax[1].set_title('Every feature: drop robustly below 0')
f.save(fig, "spread")
