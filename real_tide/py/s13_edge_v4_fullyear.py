#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""s13：v4林缘类型 × 全年日期 × 尖鼻咀潮位 的面板分析。

设计：Y_it = α_i + λ_t + Σ_k θ_k (Tide_t × Class_ik) + Σ_k Σ_m (Class_ik × Month_m) + ε_it
  α_i 点固定效应；λ_t 日期固定效应（吸收全部日期层面干扰与全区共同季节）；
  Class×Month 吸收各林缘类型不同的物候（s10显示季节振幅不同）；
  θ_k 为林缘类型k相对基准类型的潮位斜率差（基准=核心区；另以向陆林缘为基准直接检验“向海−向陆”）。
推断：日期×空间区块双向聚类；关键项再做wild cluster bootstrap（WCR，Webb权重，日期/双向）。
暴露：尖鼻咀TBT天文潮（相对平均海平面，m）；另做阈值模型 High_t = 1[Tide_t >= --high]。
输出：tidedata/s13_edge_v4_fullyear/ 下 results.csv、RESULTS_SUMMARY.md
"""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd

from s09_edge_datefe_bootstrap import project, Clusters, cov_2way, wild_p, t_p, bh_fdr, normal_p

CLASS_NAMES = {0: "core", 1: "transition", 2: "seaward", 3: "landward"}


def prepare(spectra, dates, subset):
    d = pd.read_csv(spectra)
    q = pd.read_csv(dates)
    q["date"] = pd.to_datetime(q.date).dt.strftime("%Y-%m-%d")
    if subset == "strict":
        q = q[q.qualified_strict.astype(str).str.lower() == "true"]
    d["date"] = pd.to_datetime(d.date).dt.strftime("%Y-%m-%d")
    d = d.merge(q[["date", "tide_height_m", "month"]].rename(columns={"tide_height_m": "tide"}), on="date")
    d = d[d.edge_class.isin([0, 1, 2, 3])].copy()
    d["NIR_B8"] = d.B8 / 1e4
    d["SWIR_B11"] = d.B11 / 1e4
    for k, n in CLASS_NAMES.items():
        d["is_" + n] = (d.edge_class == k).astype(float)
    return d


def add_exposure(d, expo, high):
    d = d.copy()
    x = d.tide if expo == "tide" else (d.tide >= high).astype(float)
    for n in CLASS_NAMES.values():
        d["x_" + n] = x * d["is_" + n]
    # 林缘类型×月份（以1月为基准月、核心区为基准类型；与点固定效应、日期固定效应共同识别）
    ctrl = []
    for n in ["transition", "seaward", "landward"]:
        for m in range(2, 13):
            c = f"{n}_m{m}"
            d[c] = d["is_" + n] * (d.month == m)
            ctrl.append(c)
    return d, ctrl


def fit(d, y, tests, controls, B, seed, boot_terms):
    use, yy, X, k = project(d, y, tests, controls, date_fe=True)
    n = len(yy); cl = Clusters(use)
    bread = np.linalg.inv(X.T @ X)
    b = bread @ X.T @ yy; e = yy - X @ b
    V = cov_2way(X, e, cl, bread, n, k)
    rows = []
    for j, name in enumerate(tests):
        se = float(np.sqrt(max(V[j, j], 0))); z = b[j] / se
        r = {"term": name, "estimate": b[j], "se": se, "ci95_low": b[j] - 1.96 * se, "ci95_high": b[j] + 1.96 * se,
             "p_normal": normal_p(z), "p_t_Gmin": t_p(z, min(cl.gd, cl.gb) - 1),
             "n_obs": n, "n_dates": cl.gd, "n_blocks": cl.gb}
        if B and name in boot_terms:
            for s_i, sch in enumerate(["date", "twoway"]):
                r["p_wcr_" + sch] = wild_p(yy, X, j, cl, bread, n, k, sch, B, np.random.default_rng(seed + s_i))[1]
        rows.append(r)
    return rows


def main():
    root = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--spectra", type=Path, required=True)
    ap.add_argument("--dates", type=Path, required=True, help="含 date、tide_height_m、month、qualified_strict 的合格日期表")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--B", type=int, default=1999)
    ap.add_argument("--high", type=float, default=None, help="可选高潮阈值（m，须与输入潮位基准一致）")
    ap.add_argument("--seed", type=int, default=20261017)
    args = ap.parse_args(); args.out_dir.mkdir(parents=True, exist_ok=True)

    indices = ["NDVI", "MNDWI", "EVI", "NDRE", "NIR_B8", "SWIR_B11"]
    rows = []
    specs = [  # (名称, 日期集, 点筛选, 暴露, 基准类型)
        ("main_all", "relaxed", "all", "tide", "core"),
        ("main_strict", "strict", "all", "tide", "core"),
        ("gmw_only", "relaxed", "gmw", "tide", "core"),
        ("vs_landward", "relaxed", "all", "tide", "landward"),
    ]
    if args.high is not None:
        specs.append(("high_tide_threshold", "relaxed", "all", "high", "core"))
    else:
        print("未设置同基准高潮阈值，仅跳过阈值模型；连续潮位模型保持原样。")
    cache = {}
    for name, subset, pts, expo, ref in specs:
        if subset not in cache:
            cache[subset] = prepare(args.spectra, args.dates, subset)
        d = cache[subset]
        if pts == "gmw":
            d = d[d.gmw_edge_zone >= 0]
        d, ctrl = add_exposure(d, expo, args.high)
        tests = ["x_" + n for n in CLASS_NAMES.values() if n != ref]
        boot = {"x_seaward", "x_landward"} if ref == "core" else {"x_seaward"}
        for y in indices:
            B = args.B if y in ("NDVI", "MNDWI") else 0
            for r in fit(d, y, tests, ctrl, B, args.seed, boot):
                r.update({"model": name, "index": y, "reference": ref, "exposure": expo})
                rows.append(r)
            print(f"{name:32s} {y:8s} done", flush=True)
    res = pd.DataFrame(rows)
    res["q_fdr_t_by_model_term"] = np.nan
    for (m, t), g in res.groupby(["model", "term"]):
        res.loc[g.index, "q_fdr_t_by_model_term"] = bh_fdr(g.p_t_Gmin.to_numpy())
    cols = ["model", "index", "reference", "exposure", "term", "estimate", "se", "ci95_low", "ci95_high",
            "p_normal", "p_t_Gmin", "p_wcr_date", "p_wcr_twoway", "q_fdr_t_by_model_term", "n_obs", "n_dates", "n_blocks"]
    res = res.reindex(columns=cols)
    res.to_csv(args.out_dir / "results.csv", index=False, encoding="utf-8-sig")
    (args.out_dir / "run_meta.json").write_text(json.dumps(vars(args), default=str, ensure_ascii=False, indent=2), encoding="utf-8")
    pd.set_option("display.width", 250)
    show = res[res.term.isin(["x_seaward", "x_landward", "x_transition"])]
    print(show[["model", "index", "term", "estimate", "ci95_low", "ci95_high", "p_t_Gmin", "p_wcr_date", "p_wcr_twoway",
                "q_fdr_t_by_model_term"]].to_string(index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
