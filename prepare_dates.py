#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""将GEE场景清单按独立日期汇总，并与已有FES潮位合并。"""
from pathlib import Path
import argparse
import json
import pandas as pd


def main():
    project = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--scene-list", required=True, type=Path)
    ap.add_argument("--tides", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--min-valid", type=float, default=0.80)
    ap.add_argument("--required-tiles", type=int, default=4)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    sc = pd.read_csv(args.scene_list)
    td = pd.read_csv(args.tides)
    required = {"scene_id", "overpass_utc", "cloud_pct", "mangrove_valid_frac", "platform", "tile_id", "solar_zenith"}
    missing = required - set(sc.columns)
    if missing:
        raise ValueError(f"scene list缺少字段: {sorted(missing)}")

    sc["overpass_utc"] = pd.to_datetime(sc["overpass_utc"], utc=True)
    sc["date"] = sc["overpass_utc"].dt.strftime("%Y-%m-%d")
    for c in ["cloud_pct", "mangrove_valid_frac", "solar_zenith"]:
        sc[c] = pd.to_numeric(sc[c], errors="coerce")

    daily = sc.groupby("date", as_index=False).agg(
        overpass_utc=("overpass_utc", "min"),
        n_tiles=("tile_id", "nunique"),
        tile_ids=("tile_id", lambda x: ",".join(sorted(set(x.dropna().astype(str))))),
        platform=("platform", lambda x: "+".join(sorted(set(x.dropna().astype(str))))),
        cloud_mean=("cloud_pct", "mean"),
        valid_mean=("mangrove_valid_frac", "mean"),
        valid_min=("mangrove_valid_frac", "min"),
        solar_zenith=("solar_zenith", "mean"),
    )

    td["overpass_utc"] = pd.to_datetime(td["overpass_utc"], utc=True)
    td["date"] = td["overpass_utc"].dt.strftime("%Y-%m-%d")
    tide_daily = td.groupby("date", as_index=False).agg(tide_height_m=("tide_height_m", "mean"))
    daily = daily.merge(tide_daily, on="date", how="left", validate="one_to_one")
    daily["year"] = pd.to_datetime(daily["date"]).dt.year
    daily["month"] = pd.to_datetime(daily["date"]).dt.month
    daily["valid_mean"] = daily["valid_mean"].fillna(0)  # 把空值填成0，防止报错
    # 强制要求瓦片数必须大于等于4（或者你想要的数字），并把下限设为0.6
    daily["qualified"] = (daily["n_tiles"] >= args.required_tiles) & (daily["valid_mean"] >= args.min_valid)

    qualified = daily.loc[daily["qualified"]].copy().sort_values("date")
    if qualified["tide_height_m"].isna().any():
        raise ValueError("部分合格日期未匹配到潮位")

    corr_cols = ["tide_height_m", "year", "month", "cloud_mean", "valid_mean", "valid_min", "solar_zenith"]
    corr = qualified[corr_cols].corr(method="spearman")["tide_height_m"].rename("spearman_rho").reset_index()
    corr = corr.rename(columns={"index": "variable"})

    daily.to_csv(args.out_dir / "dates_all_quality.csv", index=False, encoding="utf-8-sig")
    qualified.to_csv(args.out_dir / "analysis_dates_all.csv", index=False, encoding="utf-8-sig")
    corr.to_csv(args.out_dir / "tide_confounding_spearman.csv", index=False, encoding="utf-8-sig")

    report = {
        "scene_rows": int(len(sc)),
        "independent_dates": int(len(daily)),
        "qualified_dates": int(len(qualified)),
        "criteria": {"required_tiles": args.required_tiles, "min_valid_mean": args.min_valid},
        "tide_min_m": float(qualified["tide_height_m"].min()),
        "tide_max_m": float(qualified["tide_height_m"].max()),
        "tide_span_m": float(qualified["tide_height_m"].max() - qualified["tide_height_m"].min()),
        "dates_by_year": {str(k): int(v) for k, v in qualified.groupby("year").size().items()},
    }
    (args.out_dir / "date_preparation_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))
    print("\n潮位混叠诊断（Spearman）：")
    print(corr.to_string(index=False))


if __name__ == "__main__":
    main()
