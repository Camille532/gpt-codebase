#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""s11：水边线法估计红树林前缘潮滩高程，并判断过境时海水是否到达红树林前缘。

输入：s11_front_waterline_gee.js 导出的 front_points_daily.csv、front_bins_daily.csv，
      以及含逐日FES过境潮位的日期表（hhb_dates_all_quality.csv）。
原理：某像元在潮位高于其地面高程时为水、低于时为干。对每个像元，在观测潮位之间寻找
      使“干湿判错数”最少的阈值，即该像元高程（与FES同一基准，相对平均海平面）。
      从未湿过的像元高程只能判定为“高于最高过境潮位”（右删失），始终为水的为“低于最低潮位”。
输出：point_elevation.csv、bin_elevation_summary.csv、date_front_wetness.csv、
      front_wet_vs_tide.png、WATERLINE_SUMMARY.md
"""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd

BIN_LABELS = {1: "0-10 m", 2: "10-20 m", 3: "20-30 m", 4: "30-60 m",
              5: "60-100 m", 6: "100-200 m", 7: "200-500 m"}


def md_table(df):
    """不依赖tabulate的Markdown表格。"""
    def fmt(v):
        return f"{v:.3f}" if isinstance(v, (float, np.floating)) and np.isfinite(v) else str(v)
    head = "| " + " | ".join(df.columns) + " |"
    sep = "|" + "---|" * len(df.columns)
    body = ["| " + " | ".join(fmt(v) for v in r) + " |" for r in df.itertuples(index=False)]
    return "\n".join([head, sep] + body)


def json_val(v):
    """JSON不支持inf/NaN：删失中位数写成文字。"""
    if isinstance(v, str): return v
    if pd.isna(v): return None
    if np.isposinf(v): return "above_max_tide"
    if np.isneginf(v): return "below_min_tide"
    return float(v)


def waterline_threshold(tide, wet):
    """返回(高程估计, 删失类型, 判错率)。阈值h：tide>h为湿。"""
    order = np.argsort(tide)
    t, w = tide[order], wet[order].astype(int)
    n = len(t)
    # 候选阈值位置k：前k个(低潮)判干、其余判湿；k=0全湿，k=n全干
    dry_wrong = np.concatenate([[0], np.cumsum(w)])                 # 前k个中实际为湿的
    wet_wrong = np.concatenate([[0], np.cumsum((1 - w)[::-1])])[::-1]  # 后n-k个中实际为干的
    err = dry_wrong + wet_wrong
    best = np.flatnonzero(err == err.min())
    k = int(np.median(best).round())  # 并列时取中间位置
    rate = float(err.min() / n)
    if k == 0:
        return float(t[0]), "below_min", rate
    if k == n:
        return float(t[-1]), "above_max", rate
    return float((t[k - 1] + t[k]) / 2), "resolved", rate


def main():
    root = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--points", type=Path, default=root / "tidedata" / "s11_front" / "front_points_daily.csv")
    ap.add_argument("--bins", type=Path, default=root / "tidedata" / "s11_front" / "front_bins_daily.csv")
    ap.add_argument("--dates", type=Path, default=root / "tidedata" / "hhb_dates" / "hhb_dates_all_quality.csv")
    ap.add_argument("--analysis-dates", type=Path, default=root / "tidedata" / "hhb_dates" / "hhb_analysis_dates_all.csv")
    ap.add_argument("--out-dir", type=Path, default=root / "tidedata" / "s11_front")
    ap.add_argument("--min-obs", type=int, default=15)
    ap.add_argument("--max-error", type=float, default=0.2, help="判错率高于此值的像元视为干湿与潮位不一致")
    args = ap.parse_args(); args.out_dir.mkdir(parents=True, exist_ok=True)

    date_table = pd.read_csv(args.dates)
    datum = str(date_table["tide_datum"].dropna().iloc[0]) if "tide_datum" in date_table and date_table["tide_datum"].notna().any() else "输入潮位基准"
    tide = date_table[["date", "tide_height_m"]].rename(columns={"tide_height_m": "tide"})
    tide["date"] = pd.to_datetime(tide.date).dt.strftime("%Y-%m-%d")
    analysis = set(pd.to_datetime(pd.read_csv(args.analysis_dates).date).dt.strftime("%Y-%m-%d"))
    tmin, tmax = tide.tide.min(), tide.tide.max()

    # ---------- 逐点高程 ----------
    pts = pd.read_csv(args.points)
    pts["date"] = pd.to_datetime(pts.date).dt.strftime("%Y-%m-%d")
    pts = pts.dropna(subset=["MNDWI"]).merge(tide, on="date", how="inner")
    pts["wet"] = pts.MNDWI > 0
    rows = []
    for pid, g in pts.groupby("pt_id"):
        if len(g) < args.min_obs:
            continue
        h, kind, rate = waterline_threshold(g.tide.to_numpy(float), g.wet.to_numpy())
        rows.append({"pt_id": pid, "front_bin": int(g.front_bin.iloc[0]), "d_mangrove": g.d_mangrove.iloc[0],
                     "shore_id": g.shore_id.iloc[0], "lon": g.lon.iloc[0], "lat": g.lat.iloc[0],
                     "n_obs": len(g), "wet_frac": g.wet.mean(), "elev_m": h, "censor": kind,
                     "error_rate": rate, "consistent": rate <= args.max_error})
    pe = pd.DataFrame(rows)
    pe.to_csv(args.out_dir / "point_elevation.csv", index=False, encoding="utf-8-sig")

    def summarize(g):
        ok = g[g.consistent]
        res = ok[ok.censor == "resolved"].elev_m
        return pd.Series({
            "n_points": len(g), "n_consistent": len(ok),
            "share_above_max_tide": (ok.censor == "above_max").mean() if len(ok) else np.nan,
            "share_below_min_tide": (ok.censor == "below_min").mean() if len(ok) else np.nan,
            "elev_median_resolved": res.median(), "elev_p25_resolved": res.quantile(.25),
            "elev_p75_resolved": res.quantile(.75),
            # 含删失的中位数：把above_max当作+∞、below_min当作-∞后取中位数，仍落在可解析区间时才有意义
            "elev_median_with_censor": np.median(np.where(ok.censor == "above_max", np.inf,
                np.where(ok.censor == "below_min", -np.inf, ok.elev_m))) if len(ok) else np.nan,
            "median_error_rate": g.error_rate.median()})
    bs = pe.groupby("front_bin").apply(summarize, include_groups=False).reset_index()
    bs.insert(1, "distance", bs.front_bin.map(BIN_LABELS))
    bs.to_csv(args.out_dir / "bin_elevation_summary.csv", index=False, encoding="utf-8-sig")

    # ---------- 逐日面积口径：前缘湿像元比例 ----------
    bins = pd.read_csv(args.bins)
    bins["date"] = pd.to_datetime(bins.date).dt.strftime("%Y-%m-%d")
    bins["wet_frac"] = bins.n_wet / bins.n_valid.where(bins.n_valid > 0)
    wide = bins.pivot_table(index="date", columns="front_bin", values="wet_frac")
    valid = bins.pivot_table(index="date", columns="front_bin", values="n_valid")
    wide.columns = [f"wet_frac_bin{c}" for c in wide.columns]
    valid.columns = [f"n_valid_bin{c}" for c in valid.columns]
    dd = tide.merge(wide.reset_index(), on="date", how="left").merge(valid.reset_index(), on="date", how="left")
    dd["analysis_date"] = dd.date.isin(analysis)
    dd = dd.sort_values("tide")
    dd.to_csv(args.out_dir / "date_front_wetness.csv", index=False, encoding="utf-8-sig")

    # 前缘0–10 m湿比例>0.5 视为“海水到达红树林前缘”
    reach = dd["wet_frac_bin1"] > 0.5 if "wet_frac_bin1" in dd else pd.Series(False, index=dd.index)
    front = bs[bs.front_bin == 1].iloc[0] if (bs.front_bin == 1).any() else None
    report = {
        "tide_range_all_dates_m": [float(tmin), float(tmax)],
        "n_points_with_elevation": int(len(pe)),
        "front_0_10m": None if front is None else {k: json_val(v) for k, v in front.items()},
        "dates_front_wet_gt_50pct_all": int(reach.sum()),
        "dates_front_wet_gt_50pct_analysis": int((reach & dd.analysis_date).sum()),
        "analysis_dates": int(dd.analysis_date.sum()),
        "lowest_tide_with_front_wet_m": float(dd.loc[reach, "tide"].min()) if reach.any() else None,
    }
    (args.out_dir / "waterline_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")

    # ---------- 图 ----------
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
        plt.rcParams["axes.unicode_minus"] = False
        fig, ax = plt.subplots(figsize=(7, 4.5))
        cmap = plt.get_cmap("viridis")
        for i, b in enumerate([1, 2, 3, 4, 5, 6, 7]):
            c = f"wet_frac_bin{b}"
            if c in dd:
                ax.plot(dd.tide, dd[c], "o-", ms=3, lw=1, color=cmap(i / 6), label=BIN_LABELS[b])
        ax.set_xlabel(f"过境潮位 (m，{datum})"); ax.set_ylabel("潮滩湿像元比例")
        ax.set_title("红树林前缘外侧潮滩：过境时湿像元比例 vs 潮位")
        ax.axhline(0.5, color="grey", lw=0.8, ls="--"); ax.legend(title="距红树林", fontsize=8)
        fig.tight_layout(); fig.savefig(args.out_dir / "front_wet_vs_tide.png", dpi=200); plt.close(fig)
    except Exception as e:  # 无图形环境时不影响表格结果
        print("[警告] 作图失败：", e)

    md = ["# s11 水边线法：红树林前缘潮滩高程与过境淹没", "",
          f"过境潮位范围（{len(tide)}日期）：{tmin:.3f} 至 {tmax:.3f} m（{datum}）。", "",
          "## 各距离档潮滩高程", "", md_table(bs), "",
          "## 海水是否到达红树林前缘（0–10 m湿像元比例>0.5）", "",
          f"- 全部{len(tide)}日期中：{report['dates_front_wet_gt_50pct_all']}个",
          f"- 严格合格日期中：{report['dates_front_wet_gt_50pct_analysis']}个",
          f"- 出现前缘湿润的最低潮位：{report['lowest_tide_with_front_wet_m']}", "",
          "解释：前缘0–10 m的潮滩高程若多数“高于最高过境潮位”，说明Sentinel-2过境时海水几乎从未到达林缘，",
          "林缘指数随潮位的变化只能来自林外更远处水体通过混合像元/邻近效应进入林缘像元，而非林下淹没。"]
    (args.out_dir / "WATERLINE_SUMMARY.md").write_text("\n".join(md), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    print(bs.to_string(index=False))


if __name__ == "__main__":
    main()
