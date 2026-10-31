#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""潮位对红树林遥感反演影响——站点流程主程序。

每个站点一个配置文件（sites/<名称>.json），所有中间文件放在该站点的工作目录：
  <work_dir>/gee/      生成的GEE脚本（复制到 code.earthengine.google.com 运行）
  <work_dir>/inputs/   从Google Drive下载的GEE导出文件（文件名保持GEE导出原名，不用改）
  <work_dir>/dates/    潮位、合格日期、GEE日期列表
  <work_dir>/results/  分析结果

子命令（按顺序，详见 README.md）：
  gee       根据配置（和已有的日期表）生成全部GEE脚本
  tides     读取UHSLC h383文件，匹配UTC过境潮位
  dates     保留已有质量筛选，更新实测潮位和GEE日期列表
  lagcheck  兼容旧命令：等同dates，不扫描、不校正实测潮位
  analyze   运行本地分析（水边线、光谱面板模型、高低潮合成实验、两期变化检测模拟）
"""
from pathlib import Path
import argparse
import json
import math
import re
import subprocess
import sys
import types
import shutil
import hashlib

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
PY = ROOT / "py"
TPL = ROOT / "gee_templates"
sys.path.insert(0, str(PY))


# ---------------------------------------------------------------- 配置
def load_cfg(path):
    c = json.loads(Path(path).read_text(encoding="utf-8"))
    w, s, e, n = c["aoi"]
    c.setdefault("center", [round((w + e) / 2, 4), round((s + n) / 2, 4)])
    zone = int(math.floor((c["center"][0] + 180) / 6) + 1)
    c.setdefault("crs", f"EPSG:{32600 + zone if c['center'][1] >= 0 else 32700 + zone}")
    c.setdefault("years", [2019, 2024])
    c.setdefault("ndvi_min", 0.45)
    c.setdefault("connect_max_m", 4000)
    c.setdefault("open_water_core_m", 250)
    c.setdefault("min_valid", 0.6)
    c.setdefault("strict_valid", 0.8)
    c.setdefault("high_tide_threshold_m", None)
    c.setdefault("bootstrap_B", 1999)
    y0, y1 = c["years"]
    c.setdefault("change_sim", {"early_end": y0 + 1, "late_start": y1 - 1})
    wd = Path(c.get("work_dir", ROOT / "sites" / c["prefix"].strip("_")))
    if not wd.is_absolute():
        wd = (Path(path).resolve().parent / wd).resolve()
    c["_wd"] = wd
    for sub in ["gee", "inputs", "dates", "results"]:
        (wd / sub).mkdir(parents=True, exist_ok=True)
    return c


def p(c, *parts):
    return c["_wd"].joinpath(*parts)


def inp(c, name):
    """GEE导出文件：inputs/<prefix><name>（兼容去掉前缀的写法）。"""
    for f in [p(c, "inputs", c["prefix"] + name), p(c, "inputs", name)]:
        if f.exists():
            return f
    return p(c, "inputs", c["prefix"] + name)


def run_py(script, *args):
    cmd = [sys.executable, str(PY / script)] + [str(a) for a in args]
    print("\n>>", " ".join(cmd[1:]))
    subprocess.run(cmd, check=True)


# ---------------------------------------------------------------- GEE脚本
def js_array(dates):
    dates = sorted(dates)
    rows = [",".join(f"'{x}'" for x in dates[i:i + 8]) for i in range(0, len(dates), 8)]
    return "[\n " + ",\n ".join(rows) + "]" if dates else "[]"


GUARD = ("if ({v}.length === 0) {{ throw new Error('{msg}'); }}\n")


def cmd_gee(c, quiet=False):
    sa = c["shore_split"]
    rep = {
        "__SITE_NAME__": c["site_name"], "__PREFIX__": c["prefix"], "__ASSET_ROOT__": c["asset_root"],
        "__AOI__": json.dumps(c["aoi"]), "__CENTER__": f"{c['center'][0]}, {c['center'][1]}", "__CRS__": c["crs"],
        "__SEED__": json.dumps(c["bay_seed"]), "__SHORE_AXIS__": sa["axis"], "__SHORE_SPLIT__": str(sa["value"]),
        "__LON0__": str(c["aoi"][0]), "__LAT0__": str(c["aoi"][1]),
        "__YEAR_START__": str(c["years"][0]), "__YEAR_END__": str(c["years"][1]),
        "__NDVI_MIN__": str(c["ndvi_min"]), "__CONNECT_MAX_M__": str(c["connect_max_m"]),
        "__OPEN_WATER_CORE_M__": str(c["open_water_core_m"]),
    }
    ad = p(c, "dates", "analysis_dates.csv")
    dates = pd.read_csv(ad).date.astype(str).tolist() if ad.exists() else []
    rep["__DATES_ARRAY__"] = js_array(dates)
    hl = p(c, "dates", "tide_composite_dates.js")
    if hl.exists():
        hl_js = hl.read_text(encoding="utf-8")
    else:
        msg = "高/低潮配对日期尚未生成：先运行 python pipeline.py dates，再运行 python pipeline.py gee"
        hl_js = "var HIGH_DATES = [];\nvar LOW_DATES = [];\n" + GUARD.format(v="LOW_DATES", msg=msg)
    rep["__DATES__S15"] = hl_js
    rep["__DATES__S16"] = f"var DATES = {js_array(dates)};\n" + re.sub(r"var HIGH_DATES = \[.*?\];\n", "", hl_js, flags=re.S)
    out = []
    for f in sorted(TPL.glob("g*.js")):
        t = f.read_text(encoding="utf-8")
        t = t.replace("__DATES__", rep["__DATES__S15"] if f.name.startswith("g08") else rep["__DATES__S16"])
        for k, v in rep.items():
            if not k.startswith("__DATES__S"):
                t = t.replace(k, v)
        left = re.findall(r"__[A-Z0-9_]+__", t)
        assert not left, f"{f.name} 占位符未替换：{set(left)}"
        dst = p(c, "gee", c["prefix"] + f.name)
        dst.write_text(t, encoding="utf-8", newline="\n")
        out.append(dst.name)
    if not quiet:
        print(f"已生成 {len(out)} 个GEE脚本 → {p(c, 'gee')}")
        print(f"  投影 {c['crs']}；合格日期 {len(dates)} 个；高/低潮配对日期 {'已写入' if hl.exists() else '未生成'}")


# ---------------------------------------------------------------- 潮位
def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def saved_inputs(c):
    """Once-only copies preserve QC membership and the original GEE date groups."""
    base = p(c, "dates", "_before_real_tide")
    for relative in ["dates/scene_tides.csv", "dates/analysis_dates.csv", "dates/analysis_dates_corr.csv",
                     "dates/all_dates_corr.csv", "dates/analysis_dates_corr_strict.csv",
                     "dates/prep/analysis_dates_all.csv", "dates/prep/dates_all_quality.csv",
                     "dates/tide_composite_dates.js", "gee/" + c["prefix"] + "g08_tide_composite_experiment.js",
                     "gee/" + c["prefix"] + "g09_change_detection_sim.js"]:
        original = p(c, relative)
        target = base / relative
        if original.exists() and not target.exists():
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(original, target)
    hashes = base / "export_hashes.json"
    if not hashes.exists():
        base.mkdir(parents=True, exist_ok=True)
        hashes.write_text(json.dumps({name: file_hash(inp(c, name)) for name in ["s15_summary.csv", "s16_daily_maps.csv"]}), encoding="utf-8")
    return base


def cmd_tides(c):
    import s14_fes_lagscan as s14
    sl = inp(c, "scene_list_fullyear.csv")
    if not sl.exists():
        raise SystemExit(f"缺少{sl}；请复用原g01导出的场景清单")
    tide = c["tide"]
    if tide.get("source") != "uhslc_file":
        raise SystemExit("配置tide.source必须为uhslc_file；请安装real_tide中的站点配置")
    station = Path(tide.get("file", "inputs/h383 (1).csv"))
    if not station.is_absolute():
        station = p(c, station)
    if not station.exists():
        raise SystemExit(f"缺少实测潮位文件：{station}")
    saved_inputs(c)
    s14.run_model(types.SimpleNamespace(scene_list=sl, station_file=station,
                  out=p(c, "dates", "scene_tides.csv"), max_gap_minutes=tide.get("max_gap_minutes", 90)))


def read_dates(path):
    d = pd.read_csv(path)
    d["date"] = pd.to_datetime(d.date).dt.strftime("%Y-%m-%d")
    if d.date.isna().any() or d.date.duplicated().any():
        raise ValueError(f"日期表必须一日一行：{path}")
    return d


def bool_column(values):
    return values.astype(str).str.lower().isin(["true", "1", "1.0"])


def array_dates(text, variable):
    hit = re.search(r"var\s+" + re.escape(variable) + r"\s*=\s*\[(.*?)\]\s*;", text, flags=re.S)
    return set(re.findall(r"['\"](\d{4}-\d{2}-\d{2})['\"]", hit.group(1))) if hit else set()


def cmd_dates(c, drop_unmatched=False):
    base = saved_inputs(c)
    analysis_src = base / "dates/analysis_dates.csv"
    quality_src = base / "dates/prep/dates_all_quality.csv"
    if not analysis_src.exists() or not quality_src.exists():
        raise SystemExit("此替换版需要你原流程已有的analysis_dates.csv和prep/dates_all_quality.csv，保留原质量筛选，不重建筛选规则。")
    t = pd.read_csv(p(c, "dates", "scene_tides.csv"))
    if "tide_source" not in t or not t.tide_source.eq("UHSLC_h383_observed").all():
        raise SystemExit("scene_tides.csv尚非实测来源，请先运行tides")
    t["date"] = pd.to_datetime(t.overpass_utc, utc=True).dt.strftime("%Y-%m-%d")
    daily = t.groupby("date").tide_height_m.agg(lambda x: x.mean() if np.isfinite(x).all() else np.nan).rename("tide_height_m").reset_index()

    def replace(path):
        d = read_dates(path)
        d = d.drop(columns=[col for col in d if col.startswith("tide_")])
        d = d.merge(daily, on="date", how="left", validate="one_to_one", sort=False)
        d["tide_source"] = "UHSLC_h383_observed"
        d["tide_datum"] = "UHSLC station zero"
        return d

    a, q = replace(analysis_src), replace(quality_src)
    missing = q.loc[q.tide_height_m.isna(), ["date"]].copy()
    missing["in_original_analysis"] = missing.date.isin(a.date)
    missing.to_csv(p(c, "dates", "real_tide_missing_dates.csv"), index=False, encoding="utf-8-sig")
    if (a.tide_height_m.isna().any() or q.tide_height_m.isna().any()) and not drop_unmatched:
        raise SystemExit(f"实测缺测：全部{len(q)}日期中{len(missing)}日，原分析{len(a)}日期中{int(a.tide_height_m.isna().sum())}日。详见dates/real_tide_missing_dates.csv。接受仅使用实测覆盖日期后，为dates命令加--drop-unmatched。")
    a = a[a.tide_height_m.notna()].copy()
    all_matched = q[q.tide_height_m.notna()].copy()
    if a.empty or all_matched.empty:
        raise SystemExit("没有可用实测日期")
    if "qualified_strict" not in a:
        raise SystemExit("原分析表缺少qualified_strict，无法保留原严格筛选")
    a["qualified_strict"] = bool_column(a.qualified_strict)
    a.to_csv(p(c, "dates", "analysis_dates.csv"), index=False, encoding="utf-8-sig")
    # Keep legacy filenames as data interfaces; no time correction is performed.
    a.to_csv(p(c, "dates", "analysis_dates_corr.csv"), index=False, encoding="utf-8-sig")
    a[a.qualified_strict].to_csv(p(c, "dates", "analysis_dates_corr_strict.csv"), index=False, encoding="utf-8-sig")
    all_matched.to_csv(p(c, "dates", "all_dates_corr.csv"), index=False, encoding="utf-8-sig")
    q.to_csv(p(c, "dates", "prep", "dates_all_quality.csv"), index=False, encoding="utf-8-sig")
    prep_src = base / "dates/prep/analysis_dates_all.csv"
    if prep_src.exists():
        pp = replace(prep_src)
        pp[pp.tide_height_m.notna()].to_csv(p(c, "dates", "prep", "analysis_dates_all.csv"), index=False, encoding="utf-8-sig")
    a[a.qualified_strict][["date"]].to_csv(p(c, "dates", "_strict_dates.csv"), index=False)
    run_py("make_tide_composite_dates.py", "--dates", p(c, "dates", "analysis_dates_corr.csv"),
           "--site", c["site_name"], "--out-js", p(c, "dates", "tide_composite_dates.js"))
    hl_path = p(c, "dates", "tide_composite_dates.js")
    hl = hl_path.read_text(encoding="utf-8").replace("相对平均海平面", "UHSLC站点零点")
    hl_path.write_text(hl, encoding="utf-8")
    old_hl = base / "dates/tide_composite_dates.js"
    old_text = old_hl.read_text(encoding="utf-8") if old_hl.exists() else ""
    old_rf = base / ("gee/" + c["prefix"] + "g09_change_detection_sim.js")
    old_rf_text = old_rf.read_text(encoding="utf-8") if old_rf.exists() else ""
    low_changed = array_dates(old_rf_text, "LOW_DATES") != array_dates(hl, "LOW_DATES")
    composite_changed = any(array_dates(old_text, v) != array_dates(hl, v) for v in ["HIGH_DATES", "LOW_DATES"])
    status = {"source": "UHSLC_h383_observed", "analysis_dates": len(a), "strict_dates": int(a.qualified_strict.sum()),
              "all_matched_dates": len(all_matched), "dropped_analysis_dates": missing.loc[missing.in_original_analysis, "date"].tolist(),
              "g08_groups_changed": composite_changed, "g09_low_training_changed": low_changed,
              "original_export_hashes": json.loads((base / "export_hashes.json").read_text()),
              "g06_missing_dates": [], "g07_missing_dates": []}
    for key, filename, target in [("g06_missing_dates", "front_points_daily.csv", set(a.date)),
                                  ("g07_missing_dates", "edge_spectra_v4_fullyear.csv", set(a.date))]:
        f = inp(c, filename)
        covered = set(pd.read_csv(f, usecols=["date"]).date.astype(str)) if f.exists() else set()
        status[key] = sorted(target - covered)
    p(c, "results", "real_tide_status.json").write_text(json.dumps(status, ensure_ascii=False, indent=2), encoding="utf-8")
    cmd_gee(c)
    print(f"实测日期准备完成：分析{len(a)}日，严格{int(a.qualified_strict.sum())}日，全部有效{len(all_matched)}日。")
    print(f"g08分组改变={composite_changed}；g09低潮训练组改变={low_changed}；详见results/real_tide_status.json")


def cmd_lagcheck(c, lag=None, drop_unmatched=False):
    if lag is not None:
        raise SystemExit("实测模式不接受--lag；不进行FES时间偏移扫描或校正")
    print("lagcheck仅作为旧命令别名：准备实测日期，不扫描时间偏移。")
    cmd_dates(c, drop_unmatched)


# ---------------------------------------------------------------- 分析
def cmd_analyze(c, steps):
    dates = p(c, "dates", "analysis_dates_corr.csv")
    if not dates.exists():
        raise SystemExit("缺少实测日期表：先运行tides和dates")
    status_file = p(c, "results", "real_tide_status.json")
    if not status_file.exists():
        raise SystemExit("缺少实测日期准备记录：先运行dates")
    status = json.loads(status_file.read_text(encoding="utf-8"))
    analysis_table = pd.read_csv(dates)
    if "tide_source" not in analysis_table or not analysis_table.tide_source.eq("UHSLC_h383_observed").all():
        raise SystemExit("分析日期并非实测来源")
    for step, filename, changed_key in [("composite", "s15_summary.csv", "g08_groups_changed"), ("change", "s16_daily_maps.csv", "g09_low_training_changed")]:
        if step in steps and status[changed_key] and file_hash(inp(c, filename)) == status["original_export_hashes"][filename]:
            print(f"跳过{step}：对应GEE日期/训练组已改变，inputs/{c['prefix']}{filename}仍是原导出文件。请在GEE重跑并替换导出后再分析。")
            steps = [x for x in steps if x != step]
    if "waterline" in steps and status["g06_missing_dates"]:
        raise SystemExit("原水边线导出缺少目标日期，请补导出g06：" + str(status["g06_missing_dates"]))
    if "edge" in steps and status["g07_missing_dates"]:
        raise SystemExit("原光谱导出缺少目标日期，请补导出g07：" + str(status["g07_missing_dates"]))
    R = p(c, "results")
    if "waterline" in steps and inp(c, "front_points_daily.csv").exists():
        run_py("s11_waterline_elevation.py", "--points", inp(c, "front_points_daily.csv"), "--bins", inp(c, "front_bins_daily.csv"),
               "--dates", p(c, "dates", "all_dates_corr.csv"), "--analysis-dates", p(c, "dates", "analysis_dates_corr_strict.csv"),
               "--out-dir", R / "1_waterline", "--min-obs", 40)
    if "edge" in steps and inp(c, "edge_spectra_v4_fullyear.csv").exists():
        run_py("s13_edge_v4_fullyear.py", "--spectra", inp(c, "edge_spectra_v4_fullyear.csv"), "--dates", dates,
               "--out-dir", R / "2_edge_effects", "--B", c["bootstrap_B"],
               *(["--high", c["high_tide_threshold_m"]] if c["high_tide_threshold_m"] is not None else []))
    if "composite" in steps and inp(c, "s15_summary.csv").exists():
        out = R / "3_tide_composite"; out.mkdir(exist_ok=True)
        (out / (c["prefix"] + "s15_summary.csv")).write_bytes(inp(c, "s15_summary.csv").read_bytes())
        run_py("s15_summarize.py", "--dir", out)
    if "change" in steps and inp(c, "s16_daily_maps.csv").exists():
        cs = c["change_sim"]
        run_py("s16_change_simulation.py", "--out-dir", R / "4_change_sim",
               "--early-end", cs["early_end"], "--late-start", cs["late_start"],
               "--sites", f"{c['site_name']}|{inp(c, 's16_daily_maps.csv')}|{dates}")
    print(f"\n结果目录：{R}")



def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["gee", "tides", "dates", "lagcheck", "analyze"])
    ap.add_argument("config", type=Path, help="站点配置 JSON，如 sites/example_shankou.json")
    ap.add_argument("--lag", type=float, default=None, help="仅为旧参数保留；实测模式会拒绝时间偏移")
    ap.add_argument("--steps", default="waterline,edge,composite,change", help="analyze：要运行的分析，逗号分隔")
    ap.add_argument("--drop-unmatched", action="store_true", help="dates：仅删除缺少实测潮位的日期，原影像质量标准保持不变")
    a = ap.parse_args()
    c = load_cfg(a.config)
    print(f"站点：{c['site_name']}（前缀 {c['prefix']}）；工作目录：{c['_wd']}")
    {"gee": lambda: cmd_gee(c), "tides": lambda: cmd_tides(c), "dates": lambda: cmd_dates(c, a.drop_unmatched),
     "lagcheck": lambda: cmd_lagcheck(c, a.lag, a.drop_unmatched), "analyze": lambda: cmd_analyze(c, a.steps.split(","))}[a.command]()


if __name__ == "__main__":
    main()
