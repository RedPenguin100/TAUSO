"""Per-cohort Spearman(burden, inhibition) for each off-target feature + the on-target control.

Off-target features sit left of 0 (more burden -> worse ASO); the on-target control (blue) sits
at ~0 -- binding strength to the intended target does not predict efficacy, so the effect is
off-target-specific, not a general binding-energy artefact."""
import numpy as np
import matplotlib.pyplot as plt
import _feat as f

f.setup_style()
df = f.load()

feats = list(f.BURDEN) + list(f.CONTROL)
data = [f.per_cohort_corr(df, n) for n in feats]
fig, ax = plt.subplots(figsize=(8.5, 4.8))
bp = ax.boxplot(data, vert=False, labels=feats, showmeans=True, patch_artist=True,
                medianprops=dict(color='#222'), flierprops=dict(marker='.', ms=2, alpha=.2))
for b, name in zip(bp['boxes'], feats):
    b.set_facecolor(f.BLUE if name in f.CONTROL else f.ACCENT); b.set_alpha(.5)
ax.axvline(0, color='#888', lw=1)
for i, (d, name) in enumerate(zip(data, feats)):
    tag = '  ← control' if name in f.CONTROL else ''
    ax.text(0.62, i + 1, f'med {np.median(d):+.2f} · {100*np.mean(d<0):.0f}% neg{tag}',
            va='center', fontsize=8, color=(f.BLUE if name in f.CONTROL else '#555'))
ax.set_xlim(-0.65, 0.95); ax.invert_yaxis()
ax.set_xlabel('per-cohort Spearman (feature vs inhibition)')
ax.set_title('Off-target burden predicts bad ASOs — on-target binding does not')
f.save(fig, "cohort_corr")
