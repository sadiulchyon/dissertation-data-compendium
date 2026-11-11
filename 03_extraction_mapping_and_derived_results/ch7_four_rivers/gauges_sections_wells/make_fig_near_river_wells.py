"""Near-river wells against the reference wells and the gauge (Chapter 6B, 3 Nov 2026).
Usage: python3 make_fig_near_river_wells.py ANNUAL_CSV OUT_PNG   (ANNUAL_CSV from near_river_wells.py)
"""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, sys.argv[3] if len(sys.argv) > 3 else '.')
import thesis_style as ts
ts.apply(7)
t = pd.read_csv(sys.argv[1], index_col=0)
panels = [('SY065', 'SW251', 'Sharighat', 3.6, (2013, 2021)), ('SY061', 'SW251', 'Sharighat', 1.2, (2013, 2021)),
          ('SY062', 'SW233A', 'Jaflong Spill', 4.0, (2017, 2026)), ('SY066', 'SW233A', 'Jaflong Spill', 1.6, (2017, 2026))]
fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.2), sharex=True)
for ax, (w, gc, gn, km, ph), lab in zip(axs.flat, panels, 'abcd'):
    ref90 = t[w].loc[1990:1999].median()
    s = t[w].dropna()
    deep = s > 1.5 * ref90
    ax.axvspan(ph[0], ph[1], color=ts.SAND, alpha=0.35, lw=0)
    ax.plot(t.index, t['SY087'], color=ts.GREY, lw=0.9, label='SY087 Kanaighat (reference)')
    ax.plot(t.index, t['SY060'], color=ts.GREY, lw=0.9, ls='--', label='SY060 Jaintiapur (reference)')
    ax.plot(s.index, s, color=ts.BROWN, lw=0.8, alpha=0.6)
    ax.scatter(s.index[~deep], s[~deep], s=9, color=ts.BROWN, zorder=3, label=f'{w}, as recorded')
    if w != 'SY065' and deep.any():
        ax.scatter(s.index[deep], s[deep], s=12, facecolor='white', edgecolor=ts.BROWN, zorder=3,
                   label='deep years (> 1.5 x 1990s median)')
        ax.scatter(s.index[deep], s[deep] / 3.28, s=10, marker='x', color=ts.ORANGE, zorder=3,
                   label='deep years / 3.28 (as if feet)')
    ax.set_ylim(17 if w in ('SY062', 'SY066') else 10, 0)
    ax.set_ylabel('Feb to Apr depth to water (m)')
    ax.grid(axis='y', color=ts.GRID, lw=0.5)
    a2 = ax.twinx(); a2.spines['right'].set_visible(True); a2.spines['top'].set_visible(False)
    gg = t[gc].dropna()
    a2.plot(gg.index, gg, color=ts.WATER, lw=1.2, label=f'{gn} gauge ({gc}), Jan to Mar 10th pct')
    lo, hi = gg.min(), gg.max()
    a2.set_ylim(lo - 0.5, hi + (hi - lo) * 1.2)
    a2.set_ylabel(f'{gn} water level (m PWD)', color=ts.WATER); a2.tick_params(axis='y', colors=ts.WATER)
    ax.set_title(f'({lab}) {w}, {km} km from the {gn} gauge', loc='left', fontsize=7)
    ax.set_xlim(1964, 2027)
    if lab == 'a':
        h1, l1 = ax.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels()
    if lab == 'b':
        hb, lb = ax.get_legend_handles_labels()
handles = h1 + [hb[-2], hb[-1]] + h2
labels = l1 + [lb[-2], lb[-1]] + l2
labels = [x.replace('SY065, as recorded', 'near-river well, as recorded').replace('Sharighat gauge (SW251), Jan to Mar 10th pct', 'river gauge, dry-season level (Jan to Mar 10th percentile; right axis)') for x in labels]
fig.text(0.01, 0.005, 'Shaded: lowering phase of the gauge (Sharighat 2013 to 2021; Jaflong Spill from 2017). Depth axes inverted.', fontsize=6, color=ts.MUTED)
fig.legend(handles, labels, loc='lower center', ncol=3, fontsize=6.3, bbox_to_anchor=(0.5, -0.01))
fig.tight_layout(rect=(0, 0.08, 1, 1))
fig.savefig(sys.argv[2], dpi=300)
