"""Strong ASOs are depleted of off-target hits (strong / mid / weak within cohort)."""
import numpy as np
import matplotlib.pyplot as plt
import _feat as f

f.setup_style()
df = f.load()

cats = ['strong', 'mid', 'weak']
feats = ['rRNA (total)', 'rRNA 18S', 'RNase-H1 off-tgt']
x = np.arange(len(cats)); w = 0.26
fig, ax = plt.subplots(figsize=(7.5, 4.2))
for i, name in enumerate(feats):
    d = f.strong_weak_depletion(df, name)
    ax.bar(x + (i - 1) * w, [d[c] for c in cats], w, label=name, color=[f.GREY, f.BLUE, f.ACCENT][i])
ax.set_xticks(x); ax.set_xticklabels(['strong\n(top 25%)', 'mid', 'weak\n(bottom 25%)'])
ax.set_ylabel('fraction carrying off-target burden')
ax.set_title('Off-target hits concentrate in the weak ASOs')
ax.legend(frameon=False, fontsize=9)
f.save(fig, "depletion")
