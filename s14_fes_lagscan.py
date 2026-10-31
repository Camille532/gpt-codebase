#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""s14：为影像过境时刻匹配潮位，可使用FES2022或本地潮位站CSV。

用法（需要 .venv_tides 环境，内含 eo-tides）：
  .venv_tides\\Scripts\\python s14_fes_lagscan.py --scene-list tidedata/shankou/sk_scene_list_fullyear.csv \
      --lon 109.70 --lat 21.46 --out tidedata/shankou/sk_scene_tides_fullyear.csv
输出：每条场景的 tide_height_m（过境时刻FES，相对平均海平面）以及 tide_lag_{L}h（过境时刻+L小时的FES潮位）。
之后用 --check 模式，把影像潮滩湿润比例（s11导出的bins日表）与各偏移潮位做相关，
相关最高的偏移L≈0 说明FES潮时可靠；L明显≠0（如后海湾的+2.75 h）说明需要校正或换数据。

本地站点模式：传入 --station-csv。影像过境时间始终使用场景清单中的UTC时间；
站点时间列若无时区标记，按 --station-timezone 解释后转换为UTC，再与过境时刻匹配。
输出的 overpass_utc 始终保持UTC。相邻站点记录之间线性插值，超过 --max-gap 的间隔不插值。
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd

LAGS = np.arange(-4.0, 4.01, 0.5)


def model(times_utc, lon, lat, tide_dir):
    from eo_tides.model import model_tides
    # netCDF4在Windows下打不开含中文的绝对路径，改用相对当前目录的路径（需在tide data目录下运行）
    r = model_tides(x=lon, y=lat, time=pd.DatetimeIndex(times_utc), model="FES2022",
                    directory=str(tide_dir.resolve()))
    return np.asarray(r["tide_height"].values, float)


def station_heights(scene_times, args):
    """将站点潮位线性插值到影像过境时刻；不外推或跨越过大的数据缺口。"""
    station = pd.read_csv(args.station_csv)
    missing = {args.station_time_col, args.station_tide_col} - set(station.columns)
    if missing:
        raise ValueError(f"站点文件缺少列: {sorted(missing)}")

    raw_time = pd.to_datetime(station[args.station_time_col], errors="coerce")
    if raw_time.dt.tz is None:
        station_time = raw_time.dt.tz_localize(
            args.station_timezone, ambiguous="NaT", nonexistent="NaT"
        ).dt.tz_convert("UTC")
    else:
        station_time = raw_time.dt.tz_convert("UTC")
    level = pd.to_numeric(station[args.station_tide_col], errors="coerce")
    valid = station_time.notna() & level.notna()
    obs = pd.DataFrame({"time": station_time[valid], "level": level[valid]})
    # 同一时刻的重复记录先取均值，避免插值出现不确定结果。
    obs = obs.groupby("time", as_index=False).level.mean().sort_values("time")
    if len(obs) < 2:
        raise ValueError("站点文件至少需要两个有效的时间/潮位观测值")

    x = obs.time.astype("int64").to_numpy(dtype=np.float64) / 1e9
    y = obs.level.to_numpy(dtype=float)
    target = pd.to_datetime(scene_times, utc=True).astype("int64").to_numpy(dtype=np.float64) / 1e9
    right = np.searchsorted(x, target, side="left")
    left = right - 1
    exact = (right < len(x)) & (x[np.minimum(right, len(x) - 1)] == target)
    left[exact] = right[exact]
    right[exact] = left[exact]
    in_range = (left >= 0) & (right < len(x))
    gap = np.full(len(target), np.inf)
    gap[in_range] = x[right[in_range]] - x[left[in_range]]
    max_gap = pd.Timedelta(args.max_gap).total_seconds()
    usable = in_range & (gap <= max_gap)
    values = np.full(len(target), np.nan)
    same = usable & (left == right)
    values[same] = y[left[same]]
    between = usable & (left != right)
    fraction = (target[between] - x[left[between]]) / gap[between]
    values[between] = y[left[between]] + fraction * (y[right[between]] - y[left[between]])
    if not usable.all():
        bad = pd.to_datetime(np.asarray(scene_times)[~usable], utc=True).astype(str).tolist()
        raise ValueError(
            f"{len(bad)}个影像过境时刻无法匹配站点潮位（超出站点时间范围或观测间隔大于"
            f"{args.max_gap}），示例: {bad[:5]}"
        )
    return values


def run_model(args):
    sc = pd.read_csv(args.scene_list)
    if "overpass_utc" not in sc:
        raise ValueError("场景清单缺少 overpass_utc 列")
    t = pd.to_datetime(sc.overpass_utc, utc=True)
    if args.station_csv:
        sc["tide_height_m"] = np.round(station_heights(t, args), 4)
        sc["tide_model"] = "local_tide_station"
        sc["tide_station_file"] = str(args.station_csv)
    else:
        times = t.dt.tz_convert(None)
        allt = pd.concat([times + pd.Timedelta(hours=float(L)) for L in LAGS], ignore_index=True)
        h = model(allt, args.lon, args.lat, args.tide_dir).reshape(len(LAGS), len(t))  # 一次调用，避免重复读取分潮文件
        for i, L in enumerate(LAGS):
            sc[f"tide_lag_{L:+.1f}h"] = np.round(h[i], 4)
        sc["tide_height_m"] = sc["tide_lag_+0.0h"]
        sc["tide_model"], sc["tide_lon"], sc["tide_lat"] = "FES2022", args.lon, args.lat
    args.out.parent.mkdir(parents=True, exist_ok=True)
    sc.to_csv(args.out, index=False, encoding="utf-8-sig")
    print(f"输出 {args.out}；过境潮位范围 {sc.tide_height_m.min():.2f} 至 {sc.tide_height_m.max():.2f} m（n={len(sc)}）")


def run_check(args):
    """影像水边线（各距离档湿润比例）与潮位的Spearman相关。"""
    tides = pd.read_csv(args.out)
    tides["date"] = pd.to_datetime(tides.overpass_utc, utc=True).dt.strftime("%Y-%m-%d")
    lagcols = [c for c in tides.columns if c.startswith("tide_lag_")]
    if not lagcols and "tide_height_m" in tides:
        lagcols = ["tide_height_m"]
    if not lagcols:
        raise ValueError("潮位文件中没有 tide_lag_* 或 tide_height_m 列")
    daily = tides.groupby("date")[lagcols].mean()
    b = pd.read_csv(args.check)
    b["wet"] = b.n_wet / b.n_valid.where(b.n_valid > 0)
    if args.dates is not None:
        keep = pd.read_csv(args.dates).date.astype(str)
        b = b[b.date.isin(set(keep))]
    rows = []
    for fb, g in b.groupby("front_bin"):
        g = g.set_index("date").join(daily, how="inner")
        rows.append({"front_bin": fb, **{c.replace("tide_lag_", ""): g.wet.rank().corr(g[c].rank()) for c in lagcols}})
    res = pd.DataFrame(rows).set_index("front_bin")
    res.loc["mean"] = res.mean()
    best = res.loc["mean"].idxmax()
    print(res.round(2).to_string())
    if best == "tide_height_m":
        print(f"\n站点/过境潮位与影像水边线的平均Spearman相关：{res.loc['mean', best]:.2f}")
    else:
        print(f"\n与影像水边线最吻合的FES时间偏移：{best}（平均Spearman {res.loc['mean'].max():.2f}；"
              f"过境时刻+0.0h为 {res.loc['mean', '+0.0h']:.2f}）")
    return res


def main():
    root = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--scene-list", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--lon", type=float, default=109.70)
    ap.add_argument("--lat", type=float, default=21.46)
    ap.add_argument("--tide-dir", type=Path, default=root / "tidedata")
    ap.add_argument("--station-csv", type=Path, help="本地潮位站CSV；提供时使用站点实测潮位而非FES2022")
    ap.add_argument("--station-time-col", default="timestamp", help="站点CSV中的时间列名")
    ap.add_argument("--station-tide-col", default="tide_height_m", help="站点CSV中的潮位列名，单位m")
    ap.add_argument("--station-timezone", default="UTC", help="仅用于解释站点CSV中无时区标记的输入时间；输出统一为UTC")
    ap.add_argument("--max-gap", default="2h", help="允许插值的最大相邻观测间隔，如 30min、2h")
    ap.add_argument("--check", type=Path, help="s11导出的 front_bins_daily CSV；给出时只做检验")
    ap.add_argument("--dates", type=Path, help="检验时只用这些日期（如合格日期表）")
    args = ap.parse_args()
    run_check(args) if args.check else run_model(args)


if __name__ == "__main__":
    main()
