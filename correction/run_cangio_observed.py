#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Cangio migration entry: discover original date filenames, write observed copies.

Default site: D:/桌面/Vietnam/sites/cangio. Put h383 (1).csv in inputs, or pass
--station-file. --list-dates is read-only. --drop-unmatched explicitly permits
restricting the sample to dates with complete observed overpass coverage.
"""
from pathlib import Path
import argparse
import json
import subprocess
import sys
import pandas as pd
from tide_data import read_date_table
from sync_gee_tide_dates import read_array

SCRIPTS = Path(__file__).resolve().parent


def discover_dates(directory):
    tables = {}
    for p in sorted(directory.rglob("*.csv")):
        header = pd.read_csv(p, nrows=0)
        if {"date", "tide_height_m"}.issubset(header) and (
                "valid_mean" in header or any(c.startswith("qualified") for c in header)):
            # These columns identify the existing QC/date tables, not lag-scan reports.
            d = read_date_table(p, allow_missing=True)
            if not d.empty:
                tables[p] = d
    if not tables:
        raise ValueError(f"{directory}中未找到含date、tide_height_m和质量统计/标志的原日期表；可用--dates-dir指定目录")
    if len({p.name.lower() for p in tables}) != len(tables):
        raise ValueError("日期目录不同子文件夹内有同名表；请用--dates-dir限定到本次实际使用的日期目录")
    return tables


def choose_analysis(tables, gee_dir, explicit=None):
    if explicit:
        path = explicit.resolve()
        matches = [p for p in tables if p.resolve() == path]
        if not matches:
            raise ValueError("--analysis-dates必须是--dates-dir内发现的原日期表")
        return matches[0]
    candidates = [p for p, d in tables.items() if "analysis" in p.name.lower() and "qualified_strict" in d]
    selected = set(read_array((gee_dir / "cg_g07_spectra.js").read_text(encoding="utf-8"), "SELECTED_DATES"))
    exact = [p for p in candidates if set(tables[p].date) == selected]
    if exact:
        candidates = exact
    # Original README step 8 uses the corrected final table. Do not rank by mtime.
    corrected = [p for p in candidates if "corr" in p.stem.lower()]
    if corrected:
        candidates = corrected
    if len(candidates) != 1:
        raise ValueError("不能唯一确定原最终分析表。请加--analysis-dates指定；候选：" + ", ".join(str(p) for p in candidates))
    return candidates[0]


def choose_waterline(tables, explicit=None):
    if explicit:
        matches = [p for p in tables if p.resolve() == explicit.resolve()]
        if len(matches) != 1:
            raise ValueError("--waterline-dates必须是--dates-dir内发现的原日期表")
        return matches[0]
    candidates = [p for p in tables if "all_quality" in p.stem.lower()]
    corrected = [p for p in candidates if "corr" in p.stem.lower()]
    candidates = corrected or candidates
    return candidates[0] if len(candidates) == 1 else None


def command(script, *args):
    return [sys.executable, str(SCRIPTS / script), *map(str, args)]


def run(cmd):
    print(subprocess.list2cmdline(cmd), flush=True)
    subprocess.run(cmd, check=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--site-dir", type=Path, default=Path("D:/桌面/Vietnam/sites/cangio"))
    ap.add_argument("--dates-dir", type=Path)
    ap.add_argument("--gee-dir", type=Path)
    ap.add_argument("--station-file", type=Path)
    ap.add_argument("--analysis-dates", type=Path)
    ap.add_argument("--waterline-dates", type=Path)
    ap.add_argument("--out-dir", type=Path)
    ap.add_argument("--drop-unmatched", action="store_true")
    ap.add_argument("--list-dates", action="store_true")
    ap.add_argument("--run-local", action="store_true", help="同时运行已有输入支持的s11、s13和s16 R1/R2；不调用GEE")
    ap.add_argument("--B", type=int, default=1999)
    ap.add_argument("--early-end", type=int, default=2020)
    ap.add_argument("--late-start", type=int, default=2023)
    args = ap.parse_args()
    dates_dir = args.dates_dir or args.site_dir / "dates"
    tables = discover_dates(dates_dir)
    for p, d in tables.items():
        print(f"{p}: {len(d)}日期，严格标志={'qualified_strict' in d}")
    if args.list_dates:
        return
    gee_dir = args.gee_dir
    if gee_dir is None:
        choices = [args.site_dir / "gee", SCRIPTS]
        gee_dir = next((p for p in choices if all((p / n).exists() for n in ["cg_g06_waterline.js", "cg_g07_spectra.js",
                    "cg_g08_tide_composite_experiment.js", "cg_g09_change_detection_sim.js"])), None)
    if gee_dir is None:
        raise ValueError("没有找到cg_g06至g09脚本；请用--gee-dir指定目录")
    analysis = choose_analysis(tables, gee_dir, args.analysis_dates)
    waterline = choose_waterline(tables, args.waterline_dates)
    print("使用原最终分析日期表：", analysis)
    if waterline is None:
        print("全部质量日期表存在歧义：保留g06原DATES；可用--waterline-dates指定后重跑。")
    inputs = args.site_dir / "inputs"
    scenes = inputs / "cg_scene_list_fullyear.csv"
    station = args.station_file
    if station is None:
        choices = sorted(inputs.glob("h383*.csv"))
        if len(choices) != 1:
            raise ValueError("请把唯一的h383 CSV放入inputs，或用--station-file指定路径")
        station = choices[0]
    if not scenes.exists() or not station.exists():
        raise FileNotFoundError(f"请确认场景清单和实测潮位路径：{scenes}；{station}")
    out = args.out_dir or args.site_dir / "observed"
    out_abs = out.resolve()
    if out_abs == dates_dir.resolve() or dates_dir.resolve() in out_abs.parents:
        raise ValueError("observed输出目录应独立于原dates目录，以保留原流程结果")
    tide_path = out / "cg_scene_tides_observed.csv"
    run(command("s14_fes_lagscan.py", "--scene-list", scenes, "--station-file", station,
                "--station-format", "uhslc", "--station-timezone", "UTC", "--station-name", "UHSLC383 Vung Tau",
                "--out", tide_path, "--allow-unmatched"))
    refresh_args = ["--tides", tide_path, "--existing-dates", *tables.keys(), "--out-dir", out / "dates"]
    if args.drop_unmatched:
        refresh_args += ["--drop-unmatched"]
    run(command("refresh_observed_dates.py", *refresh_args))
    new_analysis = out / "dates" / analysis.name
    sync_args = ["--dates", new_analysis, "--gee-dir", gee_dir, "--out-dir", out / "gee"]
    if waterline:
        sync_args += ["--waterline-dates", out / "dates" / waterline.name]
    run(command("sync_gee_tide_dates.py", *sync_args))
    cmds = []
    points, bins = inputs / "cg_front_points_daily.csv", inputs / "cg_front_bins_daily.csv"
    if waterline and points.exists() and bins.exists():
        cmds.append(command("s11_waterline_elevation.py", "--points", points, "--bins", bins, "--dates",
                            out / "dates" / waterline.name, "--analysis-dates", new_analysis, "--out-dir", out / "results/waterline"))
    spectra = inputs / "cg_edge_spectra_v4_fullyear.csv"
    if spectra.exists():
        cmds.append(command("s13_edge_v4_fullyear.py", "--spectra", spectra, "--dates", new_analysis,
                            "--out-dir", out / "results/edge_effects", "--B", args.B))
    maps = inputs / "cg_s16_daily_maps.csv"
    if maps.exists():
        cmds.append(command("s16_change_simulation.py", "--sites", f"越南芹椰|{maps}|{new_analysis}",
                            "--methods", "R1", "R2", "--out-dir", out / "results/change_rules",
                            "--early-end", args.early_end, "--late-start", args.late_start))
    out.mkdir(parents=True, exist_ok=True)
    report = {"scene_list": str(scenes), "station_file": str(station), "original_analysis_dates": str(analysis),
              "original_waterline_dates": str(waterline) if waterline else None,
              "updated_date_tables": [str(out / "dates" / p.name) for p in tables],
              "detected_local_commands": [subprocess.list2cmdline(c) for c in cmds],
              "note": "Only commands whose exact expected input filenames exist are listed. See README for manual paths. GEE exports are manual."}
    (out / "migration_inputs.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    for cmd in cmds:
        if args.run_local:
            run(cmd)
        else:
            print("可运行：", subprocess.list2cmdline(cmd))
    print(f"实测更新完成：{out}；GEE重跑判断见gee/rerun_plan.json，实际输入记录见migration_inputs.json。")


if __name__ == "__main__":
    main()
