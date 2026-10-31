#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""后海湾区域研究s08：点固定效应 + 日期/空间区块双向聚类稳健推断。

仅依赖numpy/pandas。潮位在日期层面变化，逐点斜率仅作描述，主推断来自面板模型。
"""
from pathlib import Path
import argparse
import json
import math
import numpy as np
import pandas as pd
from tide_data import read_date_table, datum_label


def normal_p(z):
    return math.erfc(abs(float(z)) / math.sqrt(2.0))


def cluster_meat(X, resid, groups):
    scores = X * resid[:, None]
    meat = np.zeros((X.shape[1], X.shape[1]))
    codes, uniques = pd.factorize(pd.Series(groups).astype(str), sort=False)
    for j in range(len(uniques)):
        s = scores[codes == j].sum(axis=0)
        meat += np.outer(s, s)
    return meat, len(uniques)


def fit_within(df, ycol, terms):
    use = df[["pt_id", "date", "block_id", ycol] + terms].dropna().copy()
    # 点内去均值等价于点固定效应；时间不变变量本身不能估计，但可与潮位交互。
    cols = [ycol] + terms
    use[cols] = use[cols] - use.groupby("pt_id")[cols].transform("mean")
    y = use[ycol].to_numpy(float)
    X = use[terms].to_numpy(float)
    keep = np.nanstd(X, axis=0) > 1e-12
    X, names = X[:, keep], [n for n, k in zip(terms, keep) if k]
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    resid = y - X @ beta
    bread = np.linalg.pinv(X.T @ X)
    md, gd = cluster_meat(X, resid, use["date"])
    mb, gb = cluster_meat(X, resid, use["block_id"])
    mi, gi = cluster_meat(X, resid, use["date"].astype(str) + "|" + use["block_id"].astype(str))
    n, k = X.shape
    def corr(g):
        return (g / (g - 1)) * ((n - 1) / (n - k)) if g > 1 and n > k else 1.0
    cov = bread @ (corr(gd)*md + corr(gb)*mb - corr(gi)*mi) @ bread
    # 非正定可能由有限簇双向修正造成；对角线截到0并在报告中留档。
    se = np.sqrt(np.maximum(np.diag(cov), 0))
    z = np.divide(beta, se, out=np.full_like(beta, np.nan), where=se > 0)
    out = pd.DataFrame({"term": names, "estimate": beta, "se_2way_cluster": se, "z": z})
    out["p_normal"] = [normal_p(v) if np.isfinite(v) else np.nan for v in z]
    out["ci95_low"] = out["estimate"] - 1.96*out["se_2way_cluster"]
    out["ci95_high"] = out["estimate"] + 1.96*out["se_2way_cluster"]
    # 距离带交互模型直接给出各带总斜率及其协方差正确的置信区间。
    if "tide" in names and "tide_x_transition" in names and "tide_x_margin" in names:
        ti, tr, ma = names.index("tide"), names.index("tide_x_transition"), names.index("tide_x_margin")
        derived=[]
        for label, weights in [
            ("tide_slope_core", {ti:1.0}),
            ("tide_slope_transition", {ti:1.0,tr:1.0}),
            ("tide_slope_margin", {ti:1.0,ma:1.0})]:
            w=np.zeros(len(names))
            for j,val in weights.items(): w[j]=val
            est=float(w@beta); ses=float(math.sqrt(max(float(w@cov@w),0)))
            zz=est/ses if ses>0 else np.nan
            derived.append({"term":label,"estimate":est,"se_2way_cluster":ses,"z":zz,
                "p_normal":normal_p(zz) if np.isfinite(zz) else np.nan,
                "ci95_low":est-1.96*ses,"ci95_high":est+1.96*ses})
        out=pd.concat([out,pd.DataFrame(derived)],ignore_index=True)
    meta = {"n_obs": n, "n_points": int(use.pt_id.nunique()), "n_dates": gd,
            "n_blocks": gb, "rank": int(np.linalg.matrix_rank(X)), "n_terms": k,
            "r2_within": float(1 - (resid @ resid)/(y @ y)) if y @ y > 0 else None}
    return out, meta


def add_model_columns(df):
    d = df.copy()
    dt = pd.to_datetime(d["date"])
    d["year_c"] = dt.dt.year - dt.dt.year.mean()
    d["month_11"] = (dt.dt.month == 11).astype(float)
    d["month_12"] = (dt.dt.month == 12).astype(float)
    d["platform_B"] = d["platform"].astype(str).str.contains("2B").astype(float)
    for c in ["cloud_mean", "solar_zenith", "canopy_ht", "optical_water_freq"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["cloud_c"] = d["cloud_mean"] - d["cloud_mean"].mean()
    d["solar_c"] = d["solar_zenith"] - d["solar_zenith"].mean()
    d["canopy_c"] = d["canopy_ht"] - d["canopy_ht"].mean()
    d["waterfreq_c"] = d["optical_water_freq"] - d["optical_water_freq"].mean()
    d["tide_x_canopy"] = d["tide"] * d["canopy_c"]
    d["tide_x_waterfreq"] = d["tide"] * d["waterfreq_c"]
    d["tide_sq"] = d["tide"] ** 2
    if "edge_zone" in d.columns:
        d["edge_transition"] = (pd.to_numeric(d["edge_zone"], errors="coerce") == 1).astype(float)
        d["edge_margin"] = (pd.to_numeric(d["edge_zone"], errors="coerce") == 2).astype(float)
        d["tide_x_transition"] = d["tide"] * d["edge_transition"]
        d["tide_x_margin"] = d["tide"] * d["edge_margin"]
    return d


def point_slopes(df, ycol, min_obs=20):
    rows=[]
    for pid,g in df.groupby("pt_id"):
        g=g.dropna(subset=[ycol,"tide"])
        if len(g)<min_obs or g.tide.std()<0.1: continue
        b=np.polyfit(g.tide.to_numpy(float),g[ycol].to_numpy(float),1)[0]
        row={"pt_id":pid,"slope":b,"n":len(g),"block_id":g.block_id.iloc[0],
             "canopy_ht":g.canopy_ht.iloc[0],"optical_water_freq":g.optical_water_freq.iloc[0],
             "shore_id":g.shore_id.iloc[0]}
        if "edge_zone" in g.columns: row["edge_zone"]=g.edge_zone.iloc[0]
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    root=Path(__file__).resolve().parent
    ap=argparse.ArgumentParser()
    ap.add_argument("--spectra",type=Path,required=True)
    ap.add_argument("--dates",type=Path,default=root/"tidedata"/"hhb_dates"/"hhb_analysis_dates_all.csv")
    ap.add_argument("--out-dir",type=Path,default=root/"tidedata"/"s08_hhb_v2")
    args=ap.parse_args(); args.out_dir.mkdir(parents=True,exist_ok=True)
    df=pd.read_csv(args.spectra); dates=read_date_table(args.dates)
    dates=dates[["date","tide_height_m","valid_mean"]].rename(columns={"tide_height_m":"tide","valid_mean":"date_valid_mean"})
    df["date"]=pd.to_datetime(df.date).dt.strftime("%Y-%m-%d")
    df=df.merge(dates,on="date",how="inner",validate="many_to_one")
    numeric=["NDVI","EVI","NDRE","MNDWI","block_id","shore_id","canopy_ht","optical_water_freq","cloud_mean","solar_zenith","tide"]
    for c in numeric: df[c]=pd.to_numeric(df[c],errors="coerce")
    # GEE v3以-9999保留冠层产品缺值的林缘；分析时恢复为NaN。
    df.loc[df["canopy_ht"] <= -9000, "canopy_ht"] = np.nan
    df=add_model_columns(df)
    base=["tide","year_c","month_11","month_12","platform_B","cloud_c","solar_c"]
    models=[]
    specs=[("NDVI_main","NDVI",base),
           ("NDVI_nonlinear","NDVI",base+["tide_sq"]),
           ("NDVI_canopy_interaction","NDVI",base+["tide_x_canopy"]),
           ("NDVI_waterfreq_interaction","NDVI",base+["tide_x_waterfreq"])]
    if "edge_zone" in df.columns:
        specs.append(("NDVI_edge_interaction","NDVI",base+["tide_x_transition","tide_x_margin"]))
        specs.append(("MNDWI_edge_interaction","MNDWI",base+["tide_x_transition","tide_x_margin"]))
    for idx in ["EVI","NDRE","MNDWI"]: specs.append((idx+"_main",idx,base))
    # 敏感性：剔除有效点数最低的2024-10-15。
    subsets={"all":df,"drop_2024_10_15":df[df.date!="2024-10-15"]}
    metas={}
    for subset_name,sub in subsets.items():
        for name,y,terms in specs:
            if subset_name!="all" and name!="NDVI_main": continue
            tab,meta=fit_within(sub,y,terms)
            tab.insert(0,"model",name); tab.insert(1,"subset",subset_name); models.append(tab)
            metas[subset_name+":"+name]=meta
    result=pd.concat(models,ignore_index=True)
    result.to_csv(args.out_dir/"panel_models_2way_cluster.csv",index=False,encoding="utf-8-sig")
    for idx in ["NDVI","EVI","NDRE","MNDWI"]:
        sl=point_slopes(df,idx); sl.to_csv(args.out_dir/f"point_slopes_{idx}.csv",index=False,encoding="utf-8-sig")
    summary={"rows":int(len(df)),"points":int(df.pt_id.nunique()),"dates":int(df.date.nunique()),
             "blocks":int(df.block_id.nunique()),
             "edge_classes":sorted(df.edge_class.dropna().unique().tolist()) if "edge_class" in df.columns else [],
             "edge_zones":sorted(pd.to_numeric(df.edge_zone,errors="coerce").dropna().unique().tolist()) if "edge_zone" in df.columns else [],
             "missing_fraction":float(1-len(df)/(df.pt_id.nunique()*df.date.nunique())),"models":metas}
    (args.out_dir/"analysis_summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({k:v for k,v in summary.items() if k!="models"},ensure_ascii=False,indent=2))
    print("\n主模型潮位项：")
    print(result[result.term=="tide"][["model","subset","estimate","se_2way_cluster","p_normal","ci95_low","ci95_high"]].to_string(index=False))

if __name__=="__main__": main()
