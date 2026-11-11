"""Baseline hydroclimate for Chapter 2 (descriptive only).
Uses the first ten complete calendar years of each BWDB record. No trends, no period comparisons.
Outputs: table2_2_baseline.csv, monthly_climatology.csv
"""
import pandas as pd, numpy as np
SRC = "src/"
cfg = {
 "Lubha":       dict(f="Lubachara_River_Data.xlsx", wl="WL_326_Lubachara", q="Q_326_Lubachara", rain="Rain_CL228_Kanaighat", gw="GWT_SY087_Kanairghat",
                     wl_name="SW326 Lubachara", rain_name="CL228 Kanaighat", gw_name="SY087 Kanairghat"),
 "Sari-Gowain": dict(f="Sari_Gowain_River_Data.xlsx", wl="WL_251_Sarighat", q="Q_251_Sarighat", rain="Rain_CL116_Lallakhal", gw="GWT_SY060_Jaintiapur",
                     wl_name="SW251 Sarighat", rain_name="CL116 Lalakhal", gw_name="SY060 Jaintiapur"),
 "Piyain":      dict(f="Piyan_River_Data.xlsx", wl="WL_233A_Jaflong_Spill", q="Q_233A_Jaflong_Spill", rain="Rain_CL48_Jaflong", gw=None,
                     wl_name="SW233A Jaflong Spill", rain_name="CL48 Jaflong", gw_name=None),
 "Jadukata":    dict(f="Jadukata_River_Data.xlsx", wl="WL_131.5_Lauergarh", q="Q_131.5_Lauergarh", rain="Rain_CL49_Lauergarh", gw="GWT_SY102_Tahirpur",
                     wl_name="SW131.5 Laurergarh", rain_name="CL49 Laurergarh", gw_name="SY102 Tahirpur"),
}
def load(f, sh):
    d = pd.read_excel(SRC + f, sheet_name=sh)
    d = d.iloc[:, :2]; d.columns = ["t", "v"]
    d["t"] = pd.to_datetime(d["t"]); d = d.dropna()
    return d
def first_decade(years_ok):
    ys = sorted(years_ok)
    return ys[0], ys[0] + 9
rows, clim = [], []
for river, c in cfg.items():
    r = load(c["f"], c["rain"]); r["y"] = r.t.dt.year; r["m"] = r.t.dt.month
    cnt = r.groupby("y").size(); full = cnt[cnt >= 0.9 * 365].index
    y0, y1 = first_decade(full)
    rb = r[(r.y >= y0) & (r.y <= y1) & r.y.isin(full)]
    ann = rb.groupby("y").v.sum()
    mon = rb.groupby(["y", "m"]).v.sum().groupby("m").mean()
    jun_oct = mon.loc[6:10].sum() / mon.sum(); mar_may = mon.loc[3:5].sum() / mon.sum()
    w = load(c["f"], c["wl"]); w = w.set_index("t").v.resample("D").mean().dropna().to_frame(); w["y"] = w.index.year; w["m"] = w.index.month
    wc = w.groupby("y").size(); wfull = wc[wc >= 300].index
    wy0, wy1 = first_decade(wfull)
    wb = w[(w.y >= wy0) & (w.y <= wy1) & w.y.isin(wfull)]
    rng = wb.groupby("y").v.agg(lambda s: s.max() - s.min())
    wmon = wb.groupby("m").v.mean()
    q = load(c["f"], c["q"]); q["y"] = q.t.dt.year; q["m"] = q.t.dt.month
    qy0 = q.y.min() + (1 if q.t.min().month > 1 else 0); qy1 = qy0 + 9
    qb = q[(q.y >= qy0) & (q.y <= qy1)]
    qdry = qb[qb.m.isin([12, 1, 2, 3])].v.median(); qwet = qb[qb.m.isin([6, 7, 8, 9])].v.median()
    qmon = qb.groupby("m").v.median()
    gwd = gww = None; gy = ""
    if c["gw"]:
        g = load(c["f"], c["gw"]); g["y"] = g.t.dt.year; g["m"] = g.t.dt.month
        gb = g[(g.y >= 1990) & (g.y <= 1999)]
        gwd = gb[gb.m.isin([3, 4])].v.median(); gww = gb[gb.m.isin([8, 9, 10])].v.median(); gy = "1990-1999"
    rows.append(dict(river=river, rain_station=c["rain_name"], rain_years=f"{y0}-{y1}", rain_n_years=len(ann),
                     rain_annual_mean_mm=round(ann.mean()), rain_annual_min_mm=round(ann.min()), rain_annual_max_mm=round(ann.max()),
                     pct_jun_oct=round(100 * jun_oct), pct_mar_may=round(100 * mar_may),
                     wl_station=c["wl_name"], wl_years=f"{wy0}-{wy1}", wl_annual_range_median_m=round(rng.median(), 1),
                     wl_annual_range_min_m=round(rng.min(), 1), wl_annual_range_max_m=round(rng.max(), 1),
                     wl_low_month=int(wmon.idxmin()), wl_high_month=int(wmon.idxmax()),
                     q_years=f"{qy0}-{qy1}", q_n=len(qb), q_dry_median_m3s=round(qdry, 1), q_monsoon_median_m3s=round(qwet), q_max_measured_m3s=round(qb.v.max()),
                     gw_station=c["gw_name"], gw_years=gy, gw_depth_marapr_m=None if gwd is None else round(gwd, 1), gw_depth_augoct_m=None if gww is None else round(gww, 1)))
    for m in range(1, 13):
        clim.append(dict(river=river, month=m, rain_mm=round(mon.get(m, np.nan), 1), q_median_m3s=qmon.get(m, np.nan),
                         wl_mean_rel_m=round(wmon.get(m) - wmon.min(), 2)))
pd.DataFrame(rows).to_csv("table2_2_baseline.csv", index=False)
pd.DataFrame(clim).to_csv("monthly_climatology.csv", index=False)
print(pd.DataFrame(rows).T.to_string())
