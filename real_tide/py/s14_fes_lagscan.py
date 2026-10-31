#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Read UHSLC hourly observations. Filename/run_model retained for pipeline compatibility.
No FES model, lag scan, mean removal or automatic datum transformation is used.
"""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd


def read_observations(path):
    d = pd.read_csv(path, header=None)
    if d.shape[1] != 5 or d.empty:
        raise ValueError("h383文件必须无表头，共5列：年、月、日、UTC小时、潮位毫米值")
    d.columns = ["year", "month", "day", "hour", "height_mm"]
    for col in d:
        d[col] = pd.to_numeric(d[col], errors="raise")
    if d[["year", "month", "day", "hour"]].isna().any().any() or not d.hour.between(0, 23).all():
        raise ValueError("观测日期有空值或小时超出0—23")
    if (d[["year", "month", "day", "hour"]] % 1 != 0).any().any():
        raise ValueError("年/月/日/小时必须为整数")
    d["time"] = pd.to_datetime(d[["year", "month", "day", "hour"]], utc=True)
    good = np.isfinite(d.height_mm) & d.height_mm.ne(-32767)
    print(f"实测记录{len(d)}条，缺测/无效{int((~good).sum())}条")
    obs = d.loc[good, ["time", "height_mm"]].copy()
    if obs.empty:
        raise ValueError("没有有效观测")
    if (obs.groupby("time").height_mm.nunique() > 1).any():
        raise ValueError("同一时刻存在冲突潮位，请先核实站点/重复记录")
    obs = obs.drop_duplicates("time").sort_values("time").reset_index(drop=True)
    obs["height_m"] = obs.height_mm / 1000.
    return obs


def match_observations(times, obs, max_gap_minutes=90.):
    if not np.isfinite(max_gap_minutes) or max_gap_minutes <= 0:
        raise ValueError("max_gap_minutes必须大于0")
    target = pd.DatetimeIndex(pd.to_datetime(times, utc=True, errors="raise"))
    if target.hasnans:
        raise ValueError("场景过境时间有空值")
    origin = obs.time.iloc[0]
    x = (pd.DatetimeIndex(obs.time) - origin).total_seconds().to_numpy()
    t = (target - origin).total_seconds().to_numpy()
    pos = np.searchsorted(x, t)
    right = np.minimum(pos, len(x) - 1)
    left = np.maximum(pos - 1, 0)
    exact = (pos < len(x)) & (x[right] == t)
    inside = (pos > 0) & (pos < len(x))
    gap = x[right] - x[left]
    usable = exact | (inside & (gap <= max_gap_minutes * 60))
    h = np.interp(t, x, obs.height_m.to_numpy(float))
    h[~usable] = np.nan
    status = np.where(exact, "exact", np.where(inside, "interpolated", "outside_range"))
    status[inside & ~exact & ~usable] = "gap_too_large"
    left[exact] = right[exact]
    return pd.DataFrame({"tide_height_m": np.round(h, 6), "tide_match_status": status,
                         "tide_left_utc": pd.DatetimeIndex(obs.time).take(left),
                         "tide_right_utc": pd.DatetimeIndex(obs.time).take(right)})


def run_model(args):
    sc = pd.read_csv(args.scene_list)
    if sc.empty or "overpass_utc" not in sc:
        raise ValueError("场景清单为空或缺少overpass_utc")
    sc["overpass_utc"] = pd.to_datetime(sc.overpass_utc, utc=True, errors="raise")
    sc = sc.drop(columns=[c for c in sc if c.startswith("tide_")])
    obs = read_observations(args.station_file)
    result = match_observations(sc.overpass_utc, obs, getattr(args, "max_gap_minutes", 90.))
    for col in result:
        sc[col] = result[col].to_numpy()
    sc["tide_source"] = "UHSLC_h383_observed"
    sc["tide_datum"] = "UHSLC station zero"
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    sc.to_csv(out, index=False, encoding="utf-8-sig")
    report = {"station_file": str(args.station_file), "timezone": "UTC", "input_unit": "mm", "output_unit": "m",
              "missing_value": -32767, "datum": "UHSLC station zero", "max_gap_minutes": getattr(args, "max_gap_minutes", 90.),
              "scenes": len(sc), "unmatched_scenes": int(sc.tide_height_m.isna().sum()),
              "status_counts": sc.tide_match_status.value_counts().to_dict()}
    out.with_suffix(".report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"过境潮位已写入{out}；{sc.tide_height_m.min():.3f}至{sc.tide_height_m.max():.3f}m；缺测{report['unmatched_scenes']}景")
    return sc


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scene-list", type=Path, required=True)
    ap.add_argument("--station-file", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--max-gap-minutes", type=float, default=90.)
    run_model(ap.parse_args())


if __name__ == "__main__":
    main()
