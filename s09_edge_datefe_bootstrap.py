#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""后海湾区域研究s09：距离带差异的点+日期双向固定效应与wild cluster bootstrap。

林缘相对核心区的潮位斜率差（tide×margin）在同一日期内识别，因此可加入日期固定效应，
吸收云、气溶胶、大气校正、平台、太阳角等所有日期共同干扰；核心区斜率本身被日期固定效应吸收。
推断：双向聚类（正态近似，与s08一致）+ 施加零假设的wild cluster bootstrap（WCR，Webb六点权重），
分别按日期、按空间区块、按日期×区块双向（乘积权重）重抽样。仅依赖numpy/pandas。
"""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd

from s08_analyze import add_model_columns, normal_p

WEBB = np.array([-np.sqrt(1.5), -1.0, -np.sqrt(0.5), np.sqrt(0.5), 1.0, np.sqrt(1.5)])


def load(spectra, dates_path):
    df = pd.read_csv(spectra); dates = pd.read_csv(dates_path)
    dates = dates[["date", "tide_height_m"]].rename(columns={"tide_height_m": "tide"})
    df["date"] = pd.to_datetime(df.date).dt.strftime("%Y-%m-%d")
    df = df.merge(dates, on="date", how="inner", validate="many_to_one")
    for c in ["NDVI", "EVI", "NDRE", "MNDWI", "block_id", "shore_id", "canopy_ht",
              "optical_water_freq", "cloud_mean", "solar_zenith", "tide"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df.loc[df["canopy_ht"] <= -9000, "canopy_ht"] = np.nan
    d = add_model_columns(df)
    for c in ["cloud_c", "solar_c"]:
        d["transition_x_" + c] = d["edge_transition"] * d[c]
        d["margin_x_" + c] = d["edge_margin"] * d[c]
    return d


def project(df, ycol, tests, controls, date_fe):
    """点固定效应去均值后，用FWL把检验项与y对控制项(+日期哑变量)残差化。"""
    use = df[["pt_id", "date", "block_id", ycol] + tests + controls].dropna().reset_index(drop=True)
    Z = use[controls].to_numpy(float)
    if date_fe:
        Z = np.column_stack([Z, pd.get_dummies(use["date"]).to_numpy(float)])
    W = np.column_stack([use[[ycol] + tests].to_numpy(float), Z])
    W = W - pd.DataFrame(W).groupby(use["pt_id"].to_numpy()).transform("mean").to_numpy()
    y, X, Z = W[:, 0], W[:, 1:1 + len(tests)], W[:, 1 + len(tests):]
    if Z.shape[1]:
        rank_z = np.linalg.matrix_rank(Z)
        y = y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]
        X = X - Z @ np.linalg.lstsq(Z, X, rcond=None)[0]
    else:
        rank_z = 0
    k_full = len(tests) + rank_z + use.pt_id.nunique()  # 用于小样本修正
    return use, y, X, k_full


class Clusters:
    def __init__(self, use):
        self.d, self.gd = pd.factorize(use["date"].astype(str))[0], use["date"].nunique()
        self.b, self.gb = pd.factorize(use["block_id"].astype(str))[0], use["block_id"].nunique()
        inter = use["date"].astype(str) + "|" + use["block_id"].astype(str)
        self.i, self.gi = pd.factorize(inter)[0], inter.nunique()


def meat(scores, codes, g):
    s = np.zeros((g, scores.shape[1])); np.add.at(s, codes, scores)
    return s.T @ s


def cov_2way(X, e, cl, bread, n, k):
    def c(g): return (g / (g - 1)) * ((n - 1) / (n - k))
    sc = X * e[:, None]
    m = c(cl.gd) * meat(sc, cl.d, cl.gd) + c(cl.gb) * meat(sc, cl.b, cl.gb) - c(cl.gi) * meat(sc, cl.i, cl.gi)
    return bread @ m @ bread


def cov_1way(X, e, codes, g, bread, n, k):
    corr = (g / (g - 1)) * ((n - 1) / (n - k))
    return bread @ (corr * meat(X * e[:, None], codes, g)) @ bread


def wild_p(y, X, j, cl, bread, n, k, scheme, B, rng):
    """检验第j个系数=0：施加零假设的残差，按scheme生成bootstrap t分布，返回对称双侧p。"""
    others = [c for c in range(X.shape[1]) if c != j]
    Xr = X[:, others]
    br = np.linalg.lstsq(Xr, y, rcond=None)[0] if others else np.zeros(0)
    fit_r = Xr @ br if others else np.zeros_like(y)
    e_r = y - fit_r
    A = bread @ X.T

    def tstat(yy):
        b = A @ yy; e = yy - X @ b
        if scheme == "date": V = cov_1way(X, e, cl.d, cl.gd, bread, n, k)
        elif scheme == "block": V = cov_1way(X, e, cl.b, cl.gb, bread, n, k)
        else: V = cov_2way(X, e, cl, bread, n, k)
        return b[j] / np.sqrt(max(V[j, j], 1e-300))

    t0 = tstat(y)
    tb = np.empty(B)
    for r in range(B):
        if scheme == "date": v = rng.choice(WEBB, cl.gd)[cl.d]
        elif scheme == "block": v = rng.choice(WEBB, cl.gb)[cl.b]
        else: v = rng.choice(WEBB, cl.gd)[cl.d] * rng.choice(WEBB, cl.gb)[cl.b]
        tb[r] = tstat(fit_r + v * e_r)
    return float(t0), float(np.mean(np.abs(tb) >= abs(t0)))


def fit(df, ycol, tests, controls, date_fe, B, seed, test_term="tide_x_margin"):
    use, y, X, k = project(df, ycol, tests, controls, date_fe)
    n = len(y); cl = Clusters(use)
    bread = np.linalg.inv(X.T @ X)
    beta = bread @ X.T @ y; e = y - X @ beta
    V = cov_2way(X, e, cl, bread, n, k)
    se = np.sqrt(np.maximum(np.diag(V), 0))
    rows = []
    for j, name in enumerate(tests):
        z = beta[j] / se[j] if se[j] > 0 else np.nan
        row = {"term": name, "estimate": beta[j], "se_2way_cluster": se[j], "z": z,
               "p_normal": normal_p(z), "ci95_low": beta[j] - 1.96 * se[j], "ci95_high": beta[j] + 1.96 * se[j],
               # 双向聚类以较小的簇数决定自由度（保守）。
               "p_t_Gmin": t_p(z, min(cl.gd, cl.gb) - 1)}
        if name == test_term and B > 0:
            for s_i, scheme in enumerate(["date", "block", "twoway"]):
                rng = np.random.default_rng(seed + s_i)
                _, p = wild_p(y, X, j, cl, bread, n, k, scheme, B, rng)
                row["p_wcr_" + scheme] = p
        rows.append(row)
    meta = {"n_obs": n, "n_points": int(use.pt_id.nunique()), "n_dates": cl.gd, "n_blocks": cl.gb}
    return pd.DataFrame(rows), meta


def t_p(t, dof):
    """学生t分布双侧p值（正则化不完全贝塔函数，连分式）。"""
    if not np.isfinite(t): return np.nan
    x = dof / (dof + t * t)
    return float(betainc_reg(dof / 2.0, 0.5, x))


def betainc_reg(a, b, x):
    from math import lgamma, exp, log
    if x <= 0: return 0.0
    if x >= 1: return 1.0
    lbt = lgamma(a + b) - lgamma(a) - lgamma(b) + a * log(x) + b * log(1 - x)
    if x < (a + 1) / (a + b + 2):
        return exp(lbt) * _cf(a, b, x) / a
    return 1.0 - exp(lbt) * _cf(b, a, 1 - x) / b


def _cf(a, b, x, itmax=300, eps=3e-14):
    tiny = 1e-300
    def nz(v): return v if abs(v) > tiny else tiny
    qab, qap, qam = a + b, a + 1, a - 1
    c, d = 1.0, 1 / nz(1 - qab * x / qap); h = d
    for m in range(1, itmax + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1 / nz(1 + aa * d); c = nz(1 + aa / c); h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1 / nz(1 + aa * d); c = nz(1 + aa / c)
        delta = d * c; h *= delta
        if abs(delta - 1) < eps: break
    return h


def bh_fdr(p):
    p = np.asarray(p, float); m = len(p); order = np.argsort(p)
    q = np.empty(m); prev = 1.0
    for rank, idx in reversed(list(enumerate(order, start=1))):
        prev = min(prev, p[idx] * m / rank); q[idx] = prev
    return q


def main():
    root = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--spectra", type=Path, default=root / "tidedata" / "hhb_edge_v3_inputs" / "hhb_edge_spectra_35dates_v3.csv")
    ap.add_argument("--dates", type=Path, default=root / "tidedata" / "hhb_dates" / "hhb_analysis_dates_all.csv")
    ap.add_argument("--out-dir", type=Path, default=root / "tidedata" / "s09_hhb_edge_datefe")
    ap.add_argument("--B", type=int, default=4999)
    ap.add_argument("--seed", type=int, default=20261015)
    args = ap.parse_args(); args.out_dir.mkdir(parents=True, exist_ok=True)
    df = load(args.spectra, args.dates)

    tests = ["tide_x_transition", "tide_x_margin"]
    s08_controls = ["tide", "year_c", "month_11", "month_12", "platform_B", "cloud_c", "solar_c"]
    zone_quality = ["transition_x_cloud_c", "transition_x_solar_c", "margin_x_cloud_c", "margin_x_solar_c"]
    specs = [
        # (名称, 控制项, 是否日期固定效应)
        ("A_s08_point_fe", s08_controls, False),         # 复现s08，作为对照
        ("B_point_date_fe", [], True),                   # 主模型：点+日期双向固定效应
        ("C_point_date_fe_zone_quality", zone_quality, True),  # 允许距离带对云量/太阳角有不同响应
    ]
    subsets = {"all": df, "drop_2024_10_15": df[df.date != "2024-10-15"]}
    tabs, metas = [], {}
    for sub_name, sub in subsets.items():
        for name, controls, dfe in specs:
            if sub_name != "all" and name == "A_s08_point_fe": continue
            for y in ["NDVI", "MNDWI", "EVI", "NDRE"]:
                B = args.B if (sub_name == "all" or y in ("NDVI", "MNDWI")) else 0
                tab, meta = fit(sub, y, tests, controls, dfe, B, args.seed)
                tab.insert(0, "index", y); tab.insert(0, "subset", sub_name); tab.insert(0, "model", name)
                tabs.append(tab); metas[f"{sub_name}:{name}:{y}"] = meta
                print(f"{sub_name:16s} {name:30s} {y:6s} done", flush=True)
    res = pd.concat(tabs, ignore_index=True)
    # 四个指数之间的BH-FDR（仅对林缘交互项、基于双向wild bootstrap p）。
    res["q_fdr_wcr_twoway"] = np.nan
    for (m, s), g in res[res.term == "tide_x_margin"].groupby(["model", "subset"]):
        if g["p_wcr_twoway"].notna().all() and len(g) == 4:
            res.loc[g.index, "q_fdr_wcr_twoway"] = bh_fdr(g["p_wcr_twoway"])
    res.to_csv(args.out_dir / "edge_datefe_bootstrap.csv", index=False, encoding="utf-8-sig")
    (args.out_dir / "edge_datefe_meta.json").write_text(
        json.dumps({"B": args.B, "seed": args.seed, "models": metas}, ensure_ascii=False, indent=2), encoding="utf-8")
    cols = ["model", "subset", "index", "estimate", "ci95_low", "ci95_high", "p_normal", "p_t_Gmin",
            "p_wcr_date", "p_wcr_block", "p_wcr_twoway", "q_fdr_wcr_twoway"]
    pd.set_option("display.width", 250)
    print(res[res.term == "tide_x_margin"][cols].to_string(index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
