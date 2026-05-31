"""Per-cohort paired scatter: high-rRNA vs regular mean inhibition (below diagonal = worse)."""
import matplotlib.pyplot as plt
import _feat as f

f.setup_style()
df = f.load()

d = f.high_vs_regular(df, 'rRNA (total)')
fig, ax = plt.subplots(figsize=(5.6, 5.6))
ax.scatter(d.reg, d.high, s=d.n / 4, alpha=.30, color=f.ACCENT, edgecolors='none')
ax.plot([0, 100], [0, 100], color='#888', ls='--', lw=1)
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.set_aspect('equal')
ax.set_xlabel('regular ASOs — mean inhibition (%)')
ax.set_ylabel('high-rRNA ASOs — mean inhibition (%)')
ax.set_title(f'High-rRNA vs regular, per cohort\n'
             f'{100*(d.high<d.reg).mean():.0f}% of cohorts below the line (high-rRNA worse)')
f.save(fig, "scatter")
