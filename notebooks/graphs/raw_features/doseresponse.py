"""Dose-response: within-cohort residual inhibition by rRNA / RNaseH1 burden bin (monotonic)."""
import matplotlib.pyplot as plt
import _feat as f

f.setup_style()
df = f.load()

feats = ['rRNA (total)', 'RNase-H1 off-tgt']
fig, ax = plt.subplots(1, 2, figsize=(9, 4), sharey=True)
for a, name in zip(ax, feats):
    t = f.resid_by_bin(df, name)
    a.bar(range(len(t)), t['mean'], yerr=t['se'], color=f.SEQ[:len(t)], capsize=3, edgecolor='white')
    a.set_xticks(range(len(t))); a.set_xticklabels(t.index)
    a.axhline(0, color='#999', lw=.8); a.set_title(name); a.set_xlabel('off-target burden')
    for i, m in enumerate(t['mean']):
        a.annotate(f'{m:+.1f}', (i, m), ha='center', va='bottom' if m >= 0 else 'top', fontsize=8.5)
ax[0].set_ylabel('within-cohort residual inhibition (%)')
f.save(fig, "doseresponse")
