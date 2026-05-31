"""Dumbbell: high off-target ASOs vs their regular cohort-mates (real mean inhibition + 95% CI)."""
import matplotlib.pyplot as plt
import _feat as f

f.setup_style()
df = f.load()

feats = ['rRNA (total)', 'rRNA 18S', 'RNase-H1 off-tgt', 'off-target general', 'off-target specific']
fig, ax = plt.subplots(figsize=(8.8, 4.6))
for i, name in enumerate(feats[::-1]):
    c = f.cohort_delta(df, name); reg, hi = c['reg'], c['high']; pct = 100 * (c['delta'] < 0).mean()
    ax.errorbar(reg, i, xerr=c['reg_ci'], fmt='o', ms=9.5, color=f.GREY, ecolor=f.GREY,
                elinewidth=1.4, capsize=3, zorder=2, label='regular ASOs' if i == 0 else '')
    ax.errorbar(hi, i, xerr=c['high_ci'], fmt='o', ms=9.5, color=f.ACCENT, ecolor=f.ACCENT,
                elinewidth=1.4, capsize=3, zorder=2, label='high off-target' if i == 0 else '')
    ax.annotate(f'−{reg-hi:.1f} pts · {pct:.0f}% of cohorts', ((reg + hi) / 2, i + 0.18),
                ha='center', fontsize=8.5, color='#555')
ax.set_yticks(range(len(feats))); ax.set_yticklabels(feats[::-1])
ax.set_xlabel('mean inhibition (%)')
ax.set_title('High off-target ASOs knock down ~5 points less than their cohort-mates')
ax.legend(frameon=False, loc='lower right')
fig.text(0.012, -0.02,
         "off-target burden = predicted hybridization to abundant non-target RNA "
         "(rRNA, RNase-H1, transcriptome); 'high' = top 25% within each experiment.",
         ha='left', fontsize=8, color='#666')
f.save(fig, "dumbbell")
