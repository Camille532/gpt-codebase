#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""将GEE场景清单按独立日期汇总，并与过境潮位合并（FES或实测）。已有结果建议使用refresh_observed_dates.py保留原日期筛选。"""
from pathlib import Path
import argparse
import json
import pandas as pd
from tide_data import daily_tides


def main():
    project = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--scene-list", required=True, type=Path)
    ap.add_argument("--tides", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--min-valid", type=float, default=0.80)
    ap.add_argument("--required-tiles", type=int, default=4)
    ap.add_argument("--strict-min-valid", type=float, default=0.8)
    ap.add_argument("--max-daily-span", default="30min")
    args = ap.parse_args()
    if not 0 <= args.min_valid <= args.strict_min_valid <= 1 or args.required_tiles < 1:
        raise ValueError("要求0≤min-valid≤strict-min-valid≤1，required-tiles≥1")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    sc = pd.read_csv(args.scene_list)
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

    tide_daily = daily_tides(args.tides, args.max_daily_span)
    daily = daily.merge(tide_daily, on="date", how="left", validate="one_to_one")
    daily["year"] = pd.to_datetime(daily["date"]).dt.year
    daily["month"] = pd.to_datetime(daily["date"]).dt.month
    daily["valid_mean"] = daily["valid_mean"].fillna(0)  # 把空值填成0，防止报错
    daily["qualified"] = (daily["n_tiles"] >= args.required_tiles) & (daily["valid_mean"] >= args.min_valid)

    daily["qualified_strict"] = (daily["n_tiles"] >= args.required_tiles) & (daily["valid_mean"] >= args.strict_min_valid)

    qualified = daily.loc[daily["qualified"]].copy().sort_values("date")
    if qualified.empty:
        raise ValueError("没有符合影像质量要求的日期")
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
        "criteria": {"required_tiles": args.required_tiles, "min_valid_mean": args.min_valid, "strict_min_valid_mean": args.strict_min_valid},
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
