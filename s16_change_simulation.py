#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""s16：两期变化检测模拟——过境潮位造成的虚假红树林面积变化。

输入：GEE导出的 *_s16_daily_maps.csv（逐景、按林缘类型统计的有效像元与判为红树林像元）＋日期潮位表。
步骤：
 1. 逐景面积：A_t = Σ_k（类型k判为红树林的比例）×（类型k的v4面积），k=核心/过渡/向海/向陆/混合。
    只用红树林范围内有效像元比例≥--min-cover的日期（避免云遮挡造成的偏差）。
 2. 潮位效应：A_t(%) = β·潮位 + 年份效应 + 月份效应 + 误差（HC1稳健标准误）。
    年份效应吸收真实面积增减，月份效应吸收物候；β为每米潮位造成的面积偏差。
 3. 两期配对：所有“前期(≤--early-end年) × 后期(≥--late-start年)”的日期对，表观变化 ΔA = A_后 − A_前。
    按两期潮位差分组，比较表观变化，并与真实的年际变化（年份效应）对比。
输出：tidedata/s16_change_sim/ 下 CHANGE_SIM_RESULTS.md、daily_area.csv、pairs_summary.csv、图。
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd

CLS = {0: "核心", 1: "过渡", 2: "向海林缘", 3: "向陆林缘", 4: "混合林缘"}
METHODS = {"R1": "NDVI≥0.4且MNDWI<0", "R2": "仅NDVI≥0.4", "RF": "RF分类（低潮训练）"}


def ols_hc1(y, X):
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    e = y - X @ b
    n, k = X.shape
    bread = np.linalg.pinv(X.T @ X)
    V = bread @ (X.T * e ** 2) @ X @ bread * n / (n - k)
    return b, np.sqrt(np.diag(V))


def site_analysis(maps_csv, dates_csv, site, args):
    m = pd.read_csv(maps_csv)
    ref = m[m.date == "V4_TOTAL"].set_index("grp").n_valid.astype(float)
    m = m[m.date != "V4_TOTAL"].copy()
    m["grp"] = m.grp.astype(int)
    for c in ["n_valid", "n_R1", "n_R2", "n_RF"]:
        m[c] = m[c].astype(float)
    area = {k: ref.get(k, 0) * 0.01 for k in CLS}              # ha（10 m像元）
    tot_area = sum(area.values())
    q = pd.read_csv(dates_csv)
    q["date"] = pd.to_datetime(q.date).dt.strftime("%Y-%m-%d")
    q["year"] = pd.to_datetime(q.date).dt.year
    q["month"] = pd.to_datetime(q.date).dt.month
    q = q[["date", "tide_height_m", "year", "month"]].rename(columns={"tide_height_m": "tide"})

    rows = []
    for d, g in m[m.grp.isin(list(CLS))].groupby("date"):
        g = g.set_index("grp")
        cover = g.n_valid.sum() / sum(ref.get(k, 0) for k in CLS)
        r = {"date": d, "cover": cover}
        for meth in METHODS:
            frac = {k: g.loc[k, "n_" + meth] / g.loc[k, "n_valid"] if k in g.index and g.loc[k, "n_valid"] > 0 else np.nan for k in CLS}
            r[f"A_{meth}_ha"] = sum(frac[k] * area[k] for k in CLS)
            for k in CLS:
                r[f"A_{meth}_{CLS[k]}_ha"] = frac[k] * area[k]
        rows.append(r)
    a = pd.DataFrame(rows).merge(q, on="date")
    a = a[a.cover >= args.min_cover].copy()
    # 影像质量控制（与潮位无关）：核心区离水远，规则①在核心区判为红树林的比例应接近1；
    # 低于阈值说明整景受薄雾/薄云影响（如后海湾2019-09-27），剔除并记录
    core_frac = a[f"A_R1_{CLS[0]}_ha"] / area[0]
    bad = a[core_frac < args.qc_core_r1]
    if len(bad):
        print(f"{site}: 质控剔除 {len(bad)} 景（核心区R1比例<{args.qc_core_r1}）：", ", ".join(bad.date))
    a = a[core_frac >= args.qc_core_r1].copy()
    a.insert(0, "site", site)

    # 2. 潮位效应
    eff = []
    yrs = sorted(a.year.unique()); mons = sorted(a.month.unique())
    X = np.column_stack([a.tide.to_numpy(float)] + [(a.year == y).astype(float) for y in yrs] +
                        [(a.month == mm).astype(float) for mm in mons[1:]])
    for meth in METHODS:
        for part, col in [("全林", f"A_{meth}_ha")] + [(CLS[k], f"A_{meth}_{CLS[k]}_ha") for k in CLS]:
            base = tot_area if part == "全林" else area[[k for k in CLS if CLS[k] == part][0]]
            y = 100 * a[col].to_numpy(float) / base
            ok = np.isfinite(y)
            b, se = ols_hc1(y[ok], X[ok])
            yr = b[1:1 + len(yrs)]
            trend = np.polyfit(np.array(yrs, float), yr, 1)[0]   # 年份效应的线性趋势（%/年），代表真实变化
            eff.append({"site": site, "method": METHODS[meth], "part": part, "slope_pct_per_m": b[0], "se": se[0],
                        "ci_low": b[0] - 1.96 * se[0], "ci_high": b[0] + 1.96 * se[0],
                        "slope_ha_per_m": b[0] * base / 100, "true_trend_pct_per_yr": trend, "n_dates": int(ok.sum())})
    eff = pd.DataFrame(eff)

    # 3. 两期配对
    early = a[a.year <= args.early_end]; late = a[a.year >= args.late_start]
    pr = []
    for meth in METHODS:
        col = f"A_{meth}_ha"
        E = early[["date", "tide", "month", col]].rename(columns={"date": "d1", "tide": "t1", "month": "m1", col: "A1"})
        L = late[["date", "tide", "month", col]].rename(columns={"date": "d2", "tide": "t2", "month": "m2", col: "A2"})
        P = E.merge(L, how="cross")
        P["dtide"] = P.t2 - P.t1
        P["dA_pct"] = 100 * (P.A2 - P.A1) / tot_area
        P["same_season"] = ((P.m1 - P.m2).abs() % 12).isin([0, 1, 11])
        P["method"] = METHODS[meth]
        pr.append(P)
    P = pd.concat(pr, ignore_index=True)
    bins = [-9, -1.0, -0.3, 0.3, 1.0, 9]
    labels = ["后期低≥1 m", "后期低0.3–1 m", "潮位相近(±0.3 m)", "后期高0.3–1 m", "后期高≥1 m"]
    P["dtide_bin"] = pd.cut(P.dtide, bins, labels=labels)
    summ = (P[P.same_season].groupby(["method", "dtide_bin"], observed=True)
            .agg(n_pairs=("dA_pct", "size"), dA_mean_pct=("dA_pct", "mean"),
                 dA_p10=("dA_pct", lambda x: x.quantile(.1)), dA_p90=("dA_pct", lambda x: x.quantile(.9)))
            .reset_index())
    summ.insert(0, "site", site)
    r2 = []
    for meth, g in P[P.same_season].groupby("method"):
        c = np.corrcoef(g.dtide, g.dA_pct)[0, 1]
        r2.append({"site": site, "method": meth, "r2_dA_vs_dtide": c ** 2, "n_pairs": len(g)})
    return a, eff, summ, pd.DataFrame(r2), P, tot_area


def main():
    root = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", type=Path, default=root / "tidedata" / "s16_change_sim")
    ap.add_argument("--min-cover", type=float, default=0.8)
    ap.add_argument("--qc-core-r1", type=float, default=0.8, help="核心区R1红树林比例下限，0表示不做质控")
    ap.add_argument("--early-end", type=int, default=2020)
    ap.add_argument("--late-start", type=int, default=2023)
    ap.add_argument("--sites", nargs="*", default=None, help="site|maps_csv|dates_csv（测试用）")
    args = ap.parse_args(); args.out_dir.mkdir(parents=True, exist_ok=True)
    sites = [("后海湾", root / "tidedata/s16_change_sim/hhb_s16_daily_maps.csv", root / "tidedata/hhb_dates_fullyear/hhb_analysis_dates_relaxed06.csv"),
             ("广西山口", root / "tidedata/s16_change_sim/sk_s16_daily_maps.csv", root / "tidedata/shankou/dates/sk_analysis_dates_relaxed06_corr.csv")]
    if args.sites:
        sites = [tuple(s.split("|")) for s in args.sites]
    A, E, S, R, PP = [], [], [], [], []
    for site, mp, dp in sites:
        if not Path(mp).exists():
            print("缺少", mp); continue
        a, eff, summ, r2, P, tot = site_analysis(mp, dp, site, args)
        A.append(a); E.append(eff); S.append(summ); R.append(r2); P["site"] = site; PP.append(P)
        print(f"{site}: v4红树林 {tot:.1f} ha，可用日期 {len(a)}")
    A, E, S, R, PP = (pd.concat(x, ignore_index=True) for x in (A, E, S, R, PP))
    A.to_csv(args.out_dir / "daily_area.csv", index=False, encoding="utf-8-sig")
    E.to_csv(args.out_dir / "tide_effect_on_area.csv", index=False, encoding="utf-8-sig")
    S.to_csv(args.out_dir / "pairs_summary.csv", index=False, encoding="utf-8-sig")

    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]; plt.rcParams["axes.unicode_minus"] = False
        sites_ = PP.site.unique(); meths = list(METHODS.values())
        fig, axs = plt.subplots(len(sites_), len(meths), figsize=(5 * len(meths), 4 * len(sites_)), squeeze=False)
        for i, s in enumerate(sites_):
            for j, mth in enumerate(meths):
                g = PP[(PP.site == s) & (PP.method == mth) & PP.same_season]
                ax = axs[i, j]; ax.scatter(g.dtide, g.dA_pct, s=4, alpha=.3)
                if len(g) > 2:
                    b = np.polyfit(g.dtide, g.dA_pct, 1); xs = np.linspace(g.dtide.min(), g.dtide.max(), 10)
                    ax.plot(xs, np.polyval(b, xs), color="r"); ax.set_title(f"{s}｜{mth}\n斜率 {b[0]:+.2f}%/m", fontsize=10)
                ax.axhline(0, color="grey", lw=.6); ax.set_xlabel("两期潮位差（后期−前期，m）"); ax.set_ylabel("表观面积变化（% of v4）")
        fig.tight_layout(); fig.savefig(args.out_dir / "pair_change_vs_dtide.png", dpi=140)
    except Exception as e:
        print("[警告] 作图失败", e)

    def tbl(df, cols):
        return "\n".join(["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)] +
                         ["| " + " | ".join(f"{v:+.2f}" if isinstance(v, float) else str(v) for v in r) + " |" for r in df[cols].itertuples(index=False)])
    e1 = E[E.part.isin(["全林", "向海林缘", "向陆林缘"])].copy()
    md = ["# s16：两期变化检测模拟——潮位造成的虚假红树林面积变化", "",
          f"参数：逐景红树林有效覆盖≥{args.min_cover}；质控：核心区R1红树林比例≥{args.qc_core_r1}；前期≤{args.early_end}年，后期≥{args.late_start}年；配对限同季节（月份相差≤1）。", "",
          "## 1. 每米潮位造成的面积偏差（逐景，控制年份与月份）", "",
          tbl(e1, ["site", "method", "part", "slope_pct_per_m", "ci_low", "ci_high", "slope_ha_per_m", "true_trend_pct_per_yr", "n_dates"]), "",
          "true_trend_pct_per_yr：年份效应的线性趋势，代表同一方法测得的真实年际面积变化，用来对比潮位偏差的量级。", "",
          "## 2. 两期配对：表观面积变化（% of v4面积）按两期潮位差分组", "",
          tbl(S, ["site", "method", "dtide_bin", "n_pairs", "dA_mean_pct", "dA_p10", "dA_p90"]), "",
          "## 3. 两期潮位差能解释多少表观变化（同季节配对）", "", tbl(R, ["site", "method", "r2_dA_vs_dtide", "n_pairs"])]
    (args.out_dir / "CHANGE_SIM_RESULTS.md").write_text("\n".join(md), encoding="utf-8")
    pd.set_option("display.width", 250)
    print(e1.round(3).to_string(index=False)); print(S.round(2).to_string(index=False)); print(R.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
