"""What a 'high'-burden ASO carries: per-feature burden distribution with the top-25% cutoff."""
import numpy as np
import matplotlib.pyplot as plt
import _feat as f

f.setup_style()
df = f.load()

feats = ['rRNA (total)', 'RNase-H1 off-tgt', 'off-target general', 'off-target specific']
rng = np.random.default_rng(0)
fig, ax = plt.subplots(len(feats), 1, figsize=(8.8, 6.4))
print(f"{'feature':22s}{'%nonzero':>10s}{'high cutoff':>14s}{'typ. high':>12s}{'typ. regular':>14s}")
for a, name in zip(ax, feats):
    hi, reg = f.burden_values(df, name); s = f.burden_summary(df, name)
    for vals, col, lab in [(reg, f.GREY, 'regular'), (hi, f.ACCENT, 'high')]:
        a.scatter(np.log1p(vals), rng.uniform(-.35, .35, len(vals)), s=3, color=col, alpha=.12,
                  edgecolors='none', rasterized=True, label=lab)
    a.axvline(np.log1p(s['cutoff']), color='#333', ls='--', lw=1)
    a.set_yticks([]); a.set_ylim(-.6, .6); a.set_title(name, fontsize=10.5, pad=4)
    a.annotate(f"{s['pct_nonzero']:.0f}% nonzero", (.99, .82), xycoords='axes fraction',
               ha='right', fontsize=8.5, color='#555')
    print(f"{name:22s}{s['pct_nonzero']:9.1f}%{s['cutoff']:14.3g}{s['high']:12.3g}{s['reg']:14.3g}")
h = [plt.Line2D([], [], marker='o', ls='', ms=7, color=c) for c in (f.GREY, f.ACCENT)]
fig.legend(h, ['regular', 'high off-target'], frameon=False, fontsize=9, ncol=2,
           loc='upper right', bbox_to_anchor=(0.995, 1.005))
ax[-1].set_xlabel('off-target burden  [log1p of the feature value]')
fig.suptitle('What a "high"-burden ASO carries — per-cohort top-25% vs the rest',
             x=.012, ha='left', fontweight='bold', fontsize=12, y=1.005)
f.save(fig, "burden_magnitude")
