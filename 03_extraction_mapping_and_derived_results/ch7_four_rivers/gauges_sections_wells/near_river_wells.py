"""Chapter 6B: near-river wells against the reference wells and the river gauges (3 Nov 2026).

Usage: python3 near_river_wells.py DELIVERY_XLSX REGIONAL_GW_XLSX WL_DIR OUT_DIR

DELIVERY_XLSX     data_source/bwdb_delivery_2026-11/data_request_saadi.xlsx (sheet 'Ground Water')
REGIONAL_GW_XLSX  data_source/requested_data/Weekly_Groundwater_Table_Data.xlsx (SY102, SY087, SY060)
WL_DIR            data_source/Water_Level_Data_1990_2026 (WL_SW251, WL_SW233A CSVs)

Annual value of a well = median depth to water in February to April (years with >= 4 readings).
Annual value of a gauge = 10th percentile of daily mean water level in January to March (years with >= 60 days),
the dry-season level of Section 6B.3.1; a February to April mean is contaminated by April flash floods.

Tests
 A. SY065 (clean record, gap 2000 to 2019) against the reference wells: change 1990-1999 to
    2020-2022 at SY065 minus the change at SY087 and SY060; 95% interval by resampling years.
 B. Unit-agnostic timing test for SY061, SY062, SY066: does each well's deep-reading episode
    coincide with the lowering phase of its gauge? Years are classed by the well's own record
    ('deep' if the Feb-Apr median exceeds 1.5 x the median of its 1990-1999 values) and the
    gauge's dry-season level is compared between deep and normal years. Spearman correlation
    of well depth with gauge level, as recorded and with deep-year values divided by 3.28.
 C. Feet hypothesis check: ratio of deep-episode to 1990s values by season (a unit change
    scales the whole year; a real dry-season drawdown deepens only the dry season).
"""
import sys, os, numpy as np, pandas as pd
from scipy import stats

deliv, regional, wldir, out = sys.argv[1:5]
os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(42)

g = pd.read_excel(deliv, 'Ground Water')
g = pd.DataFrame({'well': g['Well alias'], 'date': pd.to_datetime(g['Observation date']),
                  'depth': g['Depth to water (m)']})
x = pd.ExcelFile(regional)
for sh in x.sheet_names:
    d = pd.read_excel(x, sh)
    g = pd.concat([g, pd.DataFrame({'well': sh.split()[1], 'date': pd.to_datetime(d.iloc[:, 0]),
                                    'depth': d.iloc[:, 1]})])
g = g.dropna()
g = g[g.depth >= 0]                       # one negative reading at SY060 (2007), as in 6B.3.1
g['y'] = g.date.dt.year; g['m'] = g.date.dt.month

def annual_feb_apr(s):
    s = s[s.m.isin([2, 3, 4])]
    a = s.groupby('y').depth.agg(['count', 'median'])
    return a[a['count'] >= 4]['median']

wells = sorted(g.well.unique())
ann = pd.DataFrame({w: annual_feb_apr(g[g.well == w]) for w in wells})

def gauge_annual(code):
    d = pd.read_csv(os.path.join(wldir, f'WL_{code}_1990_2026.csv'))
    d['date'] = pd.to_datetime(d['Date'], format='%m/%d/%Y %H:%M', errors='coerce')
    d = d.dropna(subset=['date'])
    day = d.groupby(d.date.dt.normalize()).WL.mean()
    day = day[day.index.month.isin([1, 2, 3])]
    a = day.groupby(day.index.year).agg(['count', lambda v: np.percentile(v, 10)])
    a.columns = ['count', 'p10']
    return a[a['count'] >= 60]['p10']

gauges = {'SW251': gauge_annual('SW251'), 'SW233A': gauge_annual('SW233A')}
pairs = {'SY061': 'SW251', 'SY065': 'SW251', 'SY062': 'SW233A', 'SY066': 'SW233A'}
dist = {'SY061': 1.2, 'SY065': 3.6, 'SY062': 4.0, 'SY066': 1.6}
tab = ann.copy()
for k, v in gauges.items():
    tab[k] = v
tab.round(3).to_csv(os.path.join(out, 'annual_feb_apr_wells_gauges.csv'), index_label='year')

lines = []
P = lines.append

# ---------- Test A: SY065 against the reference wells
P('TEST A  SY065 (3.6 km from Sharighat) against the reference wells, Feb-Apr median depth (m)')
base, late = list(range(1990, 2000)), [2020, 2021, 2022]
def did(w, ref, nboot=2000):
    b = ann.loc[ann.index.isin(base), [w, ref]].dropna()
    l = ann.loc[ann.index.isin(late), [w, ref]].dropna()
    est = (l[w].mean() - b[w].mean()) - (l[ref].mean() - b[ref].mean())
    bs = []
    for _ in range(nboot):
        bb = b.iloc[rng.integers(0, len(b), len(b))]
        ll = l.iloc[rng.integers(0, len(l), len(l))]
        bs.append((ll[w].mean() - bb[w].mean()) - (ll[ref].mean() - bb[ref].mean()))
    return est, np.percentile(bs, [2.5, 97.5]), len(b), len(l), b, l
for ref in ['SY087', 'SY060']:
    est, ci, nb, nl, b, l = did('SY065', ref)
    P(f'  SY065 1990s {b["SY065"].mean():.2f} -> 2020-22 {l["SY065"].mean():.2f} (change {l["SY065"].mean()-b["SY065"].mean():+.2f});'
      f' {ref} {b[ref].mean():.2f} -> {l[ref].mean():.2f} (change {l[ref].mean()-b[ref].mean():+.2f})')
    P(f'  excess deepening at SY065 relative to {ref}: {est:+.2f} m (95% {ci[0]:+.2f} to {ci[1]:+.2f}); years {nb} and {nl}')
# baseline coupling of SY065 with the references and its own longer record
b = ann.loc[1990:1999, ['SY065', 'SY087', 'SY060']].dropna()
for ref in ['SY087', 'SY060']:
    r, p = stats.pearsonr(b['SY065'], b[ref])
    P(f'  1990-1999 correlation SY065 with {ref}: r = {r:.2f} (p = {p:.2f}, n = {len(b)}); without 1994 (Feb-Apr median 0.28 m, an outlier):'
      f' r = {stats.pearsonr(b.drop(1994, errors="ignore")["SY065"], b.drop(1994, errors="ignore")[ref])[0]:.2f}')
P('  SY065 late years with >= 4 Feb-Apr readings: ' + ', '.join(f'{y} {v:.2f}' for y, v in ann['SY065'].loc[2020:].dropna().items()))
s65 = ann['SY065'].dropna()
P(f'  SY065 own record: 1965-1979 mean {s65.loc[1965:1979].mean():.2f}, 1980-1989 {s65.loc[1980:1989].mean():.2f},'
  f' 1990-1999 {s65.loc[1990:1999].mean():.2f}, 2020-2022 {s65.loc[2020:2022].mean():.2f}')
gs = gauges['SW251']
P(f'  Sharighat SW251 dry-season level (Jan-Mar 10th percentile): 1990s {gs.loc[1990:1999].mean():.2f} m, 2020-22 {gs.loc[2020:2022].mean():.2f} m'
  f' (change {gs.loc[2020:2022].mean()-gs.loc[1990:1999].mean():+.2f} m)')
P('')

# ---------- Test C: does a unit change explain the deep episodes?
P('TEST C  Deep episodes by season: ratio of episode median to 1990s median (a change to feet gives c. 3.3 in every season)')
episodes = {'SY061': [(1977, 1984), (2010, 2019)], 'SY062': [(2008, 2016), (2019, 2021)], 'SY066': [(2005, 2009)]}
clean = {'SY061': (1985, 2009), 'SY062': (1980, 1998), 'SY066': (1996, 2004)}
seasons = {'Feb-Apr': [2, 3, 4], 'Jun-Jul': [6, 7], 'Aug-Sep': [8, 9], 'Nov-Dec': [11, 12]}
for w, eps in episodes.items():
    s = g[g.well == w]
    c0, c1 = clean[w]
    cl = s[(s.y >= c0) & (s.y <= c1)]
    for e0, e1 in eps:
        ep = s[(s.y >= e0) & (s.y <= e1)]
        txt = []
        for nm, ms in seasons.items():
            a = ep[ep.m.isin(ms)].depth.median(); c = cl[cl.m.isin(ms)].depth.median()
            txt.append(f'{nm} {a:.2f}/{c:.2f}={a/c:.1f}')
        P(f'  {w} {e0}-{e1} vs {c0}-{c1}: ' + '; '.join(txt))
# week-to-week flips by a factor near 3 inside one season
P('  Week-to-week flips (consecutive readings <= 10 days apart whose ratio lies within 2.8 to 3.8):')
for w in ['SY061', 'SY062', 'SY066', 'SY065']:
    s = g[g.well == w].sort_values('date')
    v = s.depth.values; dt = s.date.diff().dt.days.values
    rat = np.maximum(v[1:], 1e-3) / np.maximum(v[:-1], 1e-3)
    ok = (dt[1:] <= 10) & (np.minimum(v[1:], v[:-1]) >= 1.0)
    flips = ok & (((rat >= 2.8) & (rat <= 3.8)) | ((1 / rat >= 2.8) & (1 / rat <= 3.8)))
    yrs = sorted(set(s.y.values[1:][flips]))
    P(f'    {w}: {flips.sum()} flips in {len(yrs)} years {[int(y) for y in yrs]}')
P('')

# ---------- Test B: timing against the gauges, unit-agnostic
P('TEST B  Deep-reading years against the gauge (deep = Feb-Apr median > 1.5 x the well\'s 1990-1999 median)')
rows = []
for w in ['SY061', 'SY062', 'SY066', 'SY065']:
    gc = pairs[w]
    d = pd.concat([ann[w], gauges[gc]], axis=1, keys=['well', 'gauge']).dropna()
    ref90 = ann[w].loc[1990:1999].median()
    d['deep'] = d.well > 1.5 * ref90
    # SY065 has no unit problem (Test C), so it is not rescaled
    d['well_ft'] = d.well if w == 'SY065' else np.where(d.deep, d.well / 3.28, d.well)
    deep_years = list(d.index[d.deep])
    rs, ps = stats.spearmanr(d.well, d.gauge)
    rf, pf = stats.spearmanr(d.well_ft, d.gauge)
    gd = d.loc[d.deep, 'gauge'].mean() if d.deep.any() else np.nan
    gn = d.loc[~d.deep, 'gauge'].mean()
    P(f'  {w} ({dist[w]} km from {gc}): years paired {len(d)} ({d.index.min()} to {d.index.max()}); 1990s median {ref90:.2f} m')
    P(f'    deep years: {deep_years}')
    P(f'    gauge dry-season level in deep years {gd:.2f} m, in other years {gn:.2f} m')
    P(f'    Spearman well depth vs gauge level: as recorded rho = {rs:+.2f} (p = {ps:.3f});'
      f' deep years / 3.28: rho = {rf:+.2f} (p = {pf:.3f})   [river-driven deepening gives a strong negative rho]')
    # excess over the reference well, as recorded and corrected
    for refw in ['SY087', 'SY060']:
        e = pd.concat([d, ann[refw]], axis=1, join='inner').dropna()
        b0 = e.loc[1990:1999]
        if len(b0) >= 5:
            sl, ic, r, p, _ = stats.linregress(b0[refw], b0['well_ft'])
            e['resid_ft'] = e['well_ft'] - (ic + sl * e[refw])
            r2, p2 = stats.spearmanr(e['resid_ft'], e['gauge'])
            P(f'    well (deep years / 3.28) referenced to {refw} (1990s fit r = {r:.2f}, n = {len(b0)}): residual vs gauge rho = {r2:+.2f}'
              f' (p = {p2:.3f}, n = {len(e)}); mean residual 2015 on {e.loc[2015:, "resid_ft"].mean():+.2f} m (n = {len(e.loc[2015:])})')
    rows.append(dict(well=w, gauge=gc, km=dist[w], n_years=len(d), deep_years=' '.join(map(str, deep_years)),
                     gauge_deep=round(gd, 2) if gd == gd else None, gauge_other=round(gn, 2),
                     rho_recorded=round(rs, 2), p_recorded=round(ps, 3), rho_ft=round(rf, 2), p_ft=round(pf, 3)))
pd.DataFrame(rows).to_csv(os.path.join(out, 'near_river_wells_tests.csv'), index=False)
P('')
P('Gauge lowering phases (Section 6B.4.2): Sharighat steep 2017-2021 (onset 2013); Jaflong Spill ramp from 2017.')
txt = '\n'.join(lines)
open(os.path.join(out, 'near_river_wells_report.txt'), 'w').write(txt)
print(txt)
