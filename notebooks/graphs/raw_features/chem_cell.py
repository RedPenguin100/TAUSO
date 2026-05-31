"""rRNA / off-target effect by chemistry and by cell line (median per-cohort Spearman)."""
import numpy as np
import matplotlib.pyplot as plt
import _feat as f

f.setup_style()
df = f.load()

name = 'rRNA (total)'
fig, ax = plt.subplots(1, 2, figsize=(13, 4.4))

# by chemistry, three features
chems = ['cEt', "2'MOE"]; show = ['rRNA (total)', 'RNase-H1 off-tgt', 'off-target general']
xx = np.arange(len(show)); w = 0.38
for j, ch in enumerate(chems):
    sub = df[df.chem == ch]
    meds = [np.median(f.per_cohort_corr(sub, s)) for s in show]
    ax[0].barh(xx + (0.5 - j) * w, meds, w, label=ch, color=[f.ACCENT, f.BLUE][j])
ax[0].set_yticks(xx); ax[0].set_yticklabels(show); ax[0].invert_yaxis()
ax[0].axvline(0, color='#999', lw=.8)
ax[0].set_xlabel('median per-cohort Spearman'); ax[0].set_title('by chemistry'); ax[0].legend(frameon=False)

# by cell line (rRNA total), sorted by n
nc = df.groupby('Cell_line').size().sort_values(ascending=False)
cl = [c for c in nc.index if (df.Cell_line == c).sum() >= 800]
meds = [(c, np.median(f.per_cohort_corr(df[df.Cell_line == c], name)), int((df.Cell_line == c).sum())) for c in cl]
meds = [m for m in meds if np.isfinite(m[1])]
y = np.arange(len(meds))[::-1]
ax[1].barh(y, [m[1] for m in meds], color=f.ACCENT, alpha=.8)
ax[1].set_yticks(y); ax[1].set_yticklabels([f'{m[0]} (n={m[2]:,})' for m in meds], fontsize=9)
ax[1].axvline(0, color='#999', lw=.8); ax[1].set_xlabel('median per-cohort Spearman')
ax[1].set_title(f'{name} by cell line')
f.save(fig, "chem_cell")
