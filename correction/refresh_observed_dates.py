#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Update tide columns only; keep each existing date table's rows and QC flags.

Use separately for all-quality and analysis tables. No optical QC is recomputed.
Missing observation dates stop by default; --drop-unmatched explicitly drops and
records them. Outputs must be separate from inputs.
"""
from pathlib import Path
import argparse
import json
import pandas as pd
from tide_data import daily_tides, read_date_table


def refresh(existing, new_tides, drop_unmatched=False):
    old = read_date_table(existing, allow_missing=True)
    # Remove prior tide exposure, FES lags, and datum metadata. Keep quality flags.
    cols = [c for c in old if not c.startswith("tide_") and c != "tide"]
    merged = old[cols].merge(new_tides, on="date", how="left", validate="one_to_one", sort=False)
    comparison = old[["date", "tide_height_m"]].rename(columns={"tide_height_m": "old_tide_height_m"})
    comparison = comparison.merge(merged[["date", "tide_height_m"]], on="date", validate="one_to_one")
    comparison["new_minus_old_m"] = comparison.tide_height_m - comparison.old_tide_height_m
    comparison["matched"] = merged.tide_height_m.notna().to_numpy()
    bad = merged.tide_height_m.isna()
    if bad.any() and not drop_unmatched:
        raise ValueError(f"{existing}: {int(bad.sum())}日没有完整实测过境潮位：{merged.loc[bad, 'date'].tolist()}。"
                         "默认保留日期集，故停止；如接受缩小样本，请显式使用 --drop-unmatched。")
    if drop_unmatched:
        merged = merged.loc[~bad].copy()
    if merged.empty:
        raise ValueError(f"{existing}: 更新后没有可用日期")
    return merged, comparison


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tides", type=Path, required=True, help="s14输出的实测场景潮位表")
    ap.add_argument("--existing-dates", type=Path, nargs="+", required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--drop-unmatched", action="store_true")
    ap.add_argument("--max-daily-span", default="30min")
    args = ap.parse_args()
    tides = daily_tides(args.tides, args.max_daily_span)
    if "tide_source" not in tides or not tides.tide_source.eq("observed_station").all():
        raise ValueError("此工具只接受s14生成的 observed_station 潮位表")
    names = [p.name for p in args.existing_dates]
    if len(set(names)) != len(names):
        raise ValueError("多个输入日期表同名，请分目录执行")
    staged = []
    for p in args.existing_dates:
        target = args.out_dir / p.name
        if target.resolve() == p.resolve():
            raise ValueError("输出目录必须与原表分开，避免覆盖FES版本")
        table, audit = refresh(p, tides, args.drop_unmatched)
        staged.append((p, target, table, audit))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    report = []
    for p, target, table, audit in staged:
        table.to_csv(target, index=False, encoding="utf-8-sig")
        audit.to_csv(args.out_dir / (p.stem + "_tide_comparison.csv"), index=False, encoding="utf-8-sig")
        controls = [c for c in ("year", "month", "cloud_mean", "valid_mean", "valid_min", "solar_zenith") if c in table]
        if controls:
            numeric = table[["tide_height_m"] + controls].apply(pd.to_numeric, errors="coerce")
            correlation = numeric.corr(method="spearman")["tide_height_m"].rename("spearman_rho").reset_index()
            correlation.to_csv(args.out_dir / (p.stem + "_tide_confounding_spearman.csv"), index=False, encoding="utf-8-sig")
        report.append({"input": str(p), "output": str(target), "old_dates": len(audit), "new_dates": len(table),
                       "dropped_dates": audit.loc[~audit.matched, "date"].tolist(),
                       "new_tide_min_m": float(table.tide_height_m.min()),
                       "new_tide_max_m": float(table.tide_height_m.max()),
                       "dates_by_year": {str(k): int(v) for k, v in table.groupby(pd.to_datetime(table.date).dt.year).size().items()},
                       "quality_flags_preserved": [c for c in table if c.startswith("qualified")]})
    (args.out_dir / "refresh_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
