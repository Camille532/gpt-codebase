#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Scene tide matching. --station-file selects observations; otherwise retain FES.

UHSLC h383: --station-format uhslc --station-timezone UTC.
This preset means mm, missing=-32767, datum=station zero; never MSL.
No FES time-lag correction is applied to observations. --check only diagnoses
wetness correlation when the tide source is observed.
"""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from tide_data import station_table, match_station, daily_tides

LAGS = np.arange(-4.0, 4.01, 0.5)


def model(times_utc, lon, lat, tide_dir):
    from eo_tides.model import model_tides
    r = model_tides(x=lon, y=lat, time=pd.DatetimeIndex(times_utc), model="FES2022",
                    directory=str(tide_dir.resolve()))
    return np.asarray(r["tide_height"].values, float)


def run_model(args):
    if args.scene_list is None:
        raise ValueError("匹配潮位需要 --scene-list；请用原 g01 影像清单，不能用日均潮位替代")
    sc = pd.read_csv(args.scene_list)
    if sc.empty or "overpass_utc" not in sc:
        raise ValueError("影像清单为空或缺少 overpass_utc")
    t = pd.to_datetime(sc.overpass_utc, utc=True, errors="raise")
    if t.isna().any():
        raise ValueError("影像清单过境时间有空值")
    # Strip all stale FES/lag columns if an old tide CSV is used as scene input.
    sc = sc.drop(columns=[c for c in sc if c.startswith("tide_")])
    sc["overpass_utc"] = t
    if args.station_file:
        if args.station_format == "uhslc":
            args.station_format = "ymdh"
            args.station_unit = args.station_unit or "mm"
            if args.station_unit != "mm":
                raise ValueError("UHSLC预设的单位为mm；其他单位请使用ymdh格式并显式指定")
            if args.missing_value is None:
                args.missing_value = [-32767.]
            args.datum = args.datum or "UHSLC station zero"
        if not args.station_unit:
            raise ValueError("普通站点文件必须显式提供 --station-unit m/cm/mm")
        args.missing_value = args.missing_value or []
        args.datum = args.datum or "station native datum (unspecified)"
        if not np.isfinite(args.datum_offset_m):
            raise ValueError("--datum-offset-m 必须是有限数值")
        station = station_table(args)
        matched = match_station(t, station, args.max_gap)
        for c in matched:
            sc[c] = matched[c].to_numpy()
        sc["tide_height_m"] = sc.tide_height_m.round(6)
        sc["tide_source"] = sc["tide_model"] = "observed_station"
        sc["tide_datum"] = args.datum
        sc["tide_station_file"] = str(args.station_file.resolve())
        sc["tide_station_name"] = args.station_name or args.station_file.stem
        sc["tide_input_unit"] = args.station_unit
        sc["tide_input_timezone"] = args.station_timezone or "embedded timezone"
        sc["tide_datum_offset_m"] = args.datum_offset_m
        bad = sc.tide_height_m.isna()
        if bad.any() and not args.allow_unmatched:
            raise ValueError(f"{int(bad.sum())}景无法匹配实测潮位（缺测或越界），例如 {sc.loc[bad, 'overpass_utc'].head().astype(str).tolist()}。"
                             "如需导出缺测状态供筛选，使用 --allow-unmatched；不会回填FES或外推")
        report = {"source": "observed_station", "station": sc.tide_station_name.iloc[0],
                  "datum": args.datum, "unit": "m", "max_gap": args.max_gap,
                  "valid_observations": len(station), "scenes": len(sc), "unmatched_scenes": int(bad.sum()),
                  "match_status_counts": sc.tide_match_status.value_counts().to_dict(),
                  "observation_start_utc": str(station.time.min()), "observation_end_utc": str(station.time.max()),
                  "parameters": vars(args)}
    else:
        times = t.dt.tz_convert(None)
        allt = pd.concat([times + pd.Timedelta(hours=float(L)) for L in LAGS], ignore_index=True)
        h = model(allt, args.lon, args.lat, args.tide_dir).reshape(len(LAGS), len(t))
        for i, lag in enumerate(LAGS):
            sc[f"tide_lag_{lag:+.1f}h"] = np.round(h[i], 4)
        sc["tide_height_m"] = sc["tide_lag_+0.0h"]
        sc["tide_source"] = sc["tide_model"] = "FES2022"
        sc["tide_datum"] = "FES mean sea level"
        sc["tide_lon"], sc["tide_lat"] = args.lon, args.lat
        report = {"source": "FES2022", "scenes": len(sc), "parameters": vars(args)}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    sc.to_csv(args.out, index=False, encoding="utf-8-sig")
    args.out.with_suffix(".report.json").write_text(json.dumps(report, default=str, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"输出 {args.out}；有效过境潮位 {sc.tide_height_m.min():.3f} 至 {sc.tide_height_m.max():.3f} m；缺测 {sc.tide_height_m.isna().sum()}景")


def run_check(args):
    tides = pd.read_csv(args.out)
    tides["date"] = pd.to_datetime(tides.overpass_utc, utc=True).dt.strftime("%Y-%m-%d")
    observed = "tide_source" in tides and tides.tide_source.eq("observed_station").all()
    lagcols = [] if observed else [c for c in tides if c.startswith("tide_lag_")]
    if observed or not lagcols:
        lagcols = ["tide_height_m"]
        daily = daily_tides(args.out).set_index("date")[lagcols]
    else:
        daily = tides.groupby("date")[lagcols].agg(lambda x: x.mean() if x.notna().all() else np.nan)
    b = pd.read_csv(args.check)
    b["date"] = pd.to_datetime(b.date).dt.strftime("%Y-%m-%d")
    b["wet"] = b.n_wet / b.n_valid.where(b.n_valid > 0)
    if args.dates is not None:
        keep = pd.to_datetime(pd.read_csv(args.dates).date).dt.strftime("%Y-%m-%d")
        b = b[b.date.isin(set(keep))]
    rows = []
    for fb, g in b.groupby("front_bin"):
        g = g.set_index("date").join(daily, how="inner")
        row = {"front_bin": fb}
        for c in lagcols:
            pair = g[["wet", c]].dropna()
            row[c.replace("tide_lag_", "")] = pair.wet.rank().corr(pair[c].rank()) if len(pair) >= 3 else np.nan
        rows.append(row)
    if not rows:
        raise ValueError("水边线与潮位没有可比较的距离档")
    res = pd.DataFrame(rows).set_index("front_bin")
    res.loc["mean"] = res.mean()
    print(res.round(3).to_string())
    if observed:
        print("实测潮位与影像水边线的一致性诊断；不进行FES时间偏移扫描，不自动校正实测潮位。")
    elif res.loc["mean"].notna().any():
        print("最相关的输入潮位列：", res.loc["mean"].idxmax())
    if args.check_out:
        args.check_out.parent.mkdir(parents=True, exist_ok=True)
        res.to_csv(args.check_out, encoding="utf-8-sig")
    return res


def main():
    root = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scene-list", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--lon", type=float, default=109.70)
    ap.add_argument("--lat", type=float, default=21.46)
    ap.add_argument("--tide-dir", type=Path, default=root / "tidedata")
    ap.add_argument("--station-file", "--station-csv", dest="station_file", type=Path)
    ap.add_argument("--station-format", choices=["columns", "ymdh", "uhslc"], default="columns")
    ap.add_argument("--station-time-col", default="time")
    ap.add_argument("--station-tide-col", default="tide_height_m")
    ap.add_argument("--station-timezone", help="无时区时间的实际时区；h383使用UTC")
    ap.add_argument("--station-time-format")
    ap.add_argument("--station-unit", choices=["m", "cm", "mm"])
    ap.add_argument("--station-name")
    ap.add_argument("--station-sheet", help="Excel工作表名称；默认第一张表")
    ap.add_argument("--station-sep", default=",")
    ap.add_argument("--station-encoding", default="utf-8-sig")
    ap.add_argument("--missing-value", type=float, nargs="+", help="原始单位下的缺测码；UHSLC预设为-32767")
    ap.add_argument("--datum", help="输出高程基准名称；默认保留站点原始基准")
    ap.add_argument("--datum-offset-m", type=float, default=0., help="输出米值=原始值×单位换算系数+此偏移；不自动去均值")
    ap.add_argument("--max-gap", default="90min", help="允许插值的相邻有效观测间隔；默认不跨过缺失小时")
    ap.add_argument("--allow-unmatched", action="store_true")
    ap.add_argument("--check", type=Path, help="只做已有潮位表和水边线湿润比例的诊断")
    ap.add_argument("--check-out", type=Path)
    ap.add_argument("--dates", type=Path)
    args = ap.parse_args()
    run_check(args) if args.check else run_model(args)


if __name__ == "__main__":
    main()
