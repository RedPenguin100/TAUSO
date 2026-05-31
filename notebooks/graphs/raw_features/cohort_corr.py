"""Per-cohort Spearman(burden, inhibition) distribution for each off-target feature."""
import numpy as np
import matplotlib.pyplot as plt
import _feat as f

f.setup_style()
df = f.load()

feats = list(f.BURDEN)
data = [f.per_cohort_corr(df, n) for n in feats]
fig, ax = plt.subplots(figsize=(8.5, 4.6))
bp = ax.boxplot(data, vert=False, labels=feats, showmeans=True, patch_artist=True,
                medianprops=dict(color='#222'), flierprops=dict(marker='.', ms=2, alpha=.2))
for b in bp['boxes']:
    b.set_facecolor(f.ACCENT); b.set_alpha(.5)
ax.axvline(0, color='#888', lw=1)
for i, d in enumerate(data):
    ax.text(0.62, i + 1, f'med {np.median(d):+.2f} · {100*np.mean(d<0):.0f}% neg',
            va='center', fontsize=8, color='#555')
ax.set_xlim(-0.65, 0.9); ax.invert_yaxis()
ax.set_xlabel('per-cohort Spearman (burden vs inhibition)')
ax.set_title('Burden vs efficacy, per cohort')
f.save(fig, "cohort_corr")
