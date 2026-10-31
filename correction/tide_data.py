"""Shared tide IO: UTC matching, explicit units/datum, and validated date tables."""
from pathlib import Path
import numpy as np
import pandas as pd


def utc_times(values, timezone=None, time_format=None):
    # Parse separately so timezone-naive values never silently become UTC.
    result = []
    for value in values:
        stamp = pd.to_datetime(value, format=time_format, errors="raise")
        if pd.isna(stamp):
            raise ValueError("时间列有空值")
        if stamp.tzinfo is None:
            if not timezone:
                raise ValueError("无时区时间必须指定 --station-timezone，例如 Asia/Ho_Chi_Minh、Asia/Shanghai 或 UTC")
            stamp = stamp.tz_localize(timezone, ambiguous="raise", nonexistent="raise")
        result.append(stamp.tz_convert("UTC"))
    return pd.DatetimeIndex(result)


def station_table(args):
    path = args.station_file
    if args.station_format == "ymdh":
        raw = pd.read_csv(path, header=None, sep=args.station_sep, encoding=args.station_encoding)
        if raw.shape[1] != 5:
            raise ValueError("ymdh 无表头格式必须有5列：年、月、日、小时、潮位")
        raw.columns = ["year", "month", "day", "hour", "level"]
        # The uploaded h383 file uses whole UTC hours.
        for col in ["year", "month", "day", "hour"]:
            raw[col] = pd.to_numeric(raw[col], errors="raise")
            if raw[col].isna().any() or (raw[col] % 1 != 0).any():
                raise ValueError(f"{col} 必须是非空整数")
        if not raw.hour.between(0, 23).all():
            raise ValueError("ymdh 小时必须在0至23之间")
        if not args.station_timezone:
            raise ValueError("ymdh 格式必须明确指定 --station-timezone UTC（或实际时区）")
        times = pd.DatetimeIndex(pd.to_datetime(raw[["year", "month", "day", "hour"]])).tz_localize(
            args.station_timezone, ambiguous="raise", nonexistent="raise").tz_convert("UTC")
        levels = pd.to_numeric(raw.level, errors="raise").to_numpy(float)
    elif path.suffix.lower() in (".xlsx", ".xls"):
        sheet = args.station_sheet
        raw = pd.read_excel(path, sheet_name=0 if sheet is None else sheet)
    else:
        raw = pd.read_csv(path, sep=args.station_sep, encoding=args.station_encoding)
    if args.station_format != "ymdh":
        missing = {args.station_time_col, args.station_tide_col} - set(raw.columns)
        if missing:
            raise ValueError(f"实测潮位文件缺少列 {sorted(missing)}，实际列为 {raw.columns.tolist()}")
        times = utc_times(raw[args.station_time_col], args.station_timezone, args.station_time_format)
        levels = pd.to_numeric(raw[args.station_tide_col], errors="raise").to_numpy(float)
    if raw.empty:
        raise ValueError("实测潮位文件没有记录")
    # Remove only explicitly configured sentinels/empty cells; their time gaps
    # still prevent interpolation via max-gap. Never treat a sentinel as water.
    good = np.isfinite(levels) & ~np.isin(levels, args.missing_value)
    print(f"实测记录 {len(raw)} 条，缺测/无效 {int((~good).sum())} 条；这些记录不参与插值")
    times, levels = times[good], levels[good]
    if not len(levels):
        raise ValueError("实测文件没有有效潮位")
    levels = levels * {"m": 1., "cm": .01, "mm": .001}[args.station_unit] + args.datum_offset_m
    if not np.isfinite(levels).all():
        raise ValueError("潮位单位/基准换算结果非有限数值")
    station = pd.DataFrame({"time": times, "level": levels}).sort_values("time")
    duplicates = station.groupby("time")["level"].nunique()
    if (duplicates > 1).any():
        raise ValueError("同一观测时间有不同潮位，可能混入多个站点；请先选择一个代表站并处理重复值")
    return station.drop_duplicates("time").reset_index(drop=True)


def match_station(target_times, station, max_gap):
    """Linear interpolation only between close records; exact matches always valid."""
    max_gap = pd.Timedelta(max_gap)
    if pd.isna(max_gap) or max_gap <= pd.Timedelta(0):
        raise ValueError("--max-gap 必须是正时间间隔，例如 90min")
    targets = pd.DatetimeIndex(pd.to_datetime(target_times, utc=True))
    if targets.hasnans:
        raise ValueError("过境时间含空值")
    # Use relative seconds; independent of pandas datetime storage precision.
    origin = station.time.iloc[0]
    x = (pd.DatetimeIndex(station.time) - origin).total_seconds().to_numpy()
    t = (targets - origin).total_seconds().to_numpy()
    y = station.level.to_numpy(float)
    pos = np.searchsorted(x, t, side="left")
    safe = np.minimum(pos, len(x) - 1)
    exact = (pos < len(x)) & (x[safe] == t)
    left = np.maximum(pos - 1, 0)
    right = safe
    gap = x[right] - x[left]
    inside = (pos > 0) & (pos < len(x))
    valid = exact | (inside & (gap <= max_gap.total_seconds()))
    heights = np.interp(t, x, y)
    heights[~valid] = np.nan
    left[exact] = right[exact]
    status = np.where(exact, "exact", np.where(inside, "interpolated", "outside_range"))
    status[inside & ~exact & ~valid] = "gap_too_large"
    return pd.DataFrame({
        "tide_height_m": heights,
        "tide_match_status": status,
        "tide_left_utc": pd.DatetimeIndex(station.time).take(left),
        "tide_right_utc": pd.DatetimeIndex(station.time).take(right),
        "tide_bracket_minutes": np.where(exact, 0., gap / 60),
    })


def read_date_table(path, allow_missing=False):
    d = pd.read_csv(path)
    missing = {"date", "tide_height_m"} - set(d)
    if missing:
        raise ValueError(f"{path} 缺少 {sorted(missing)}")
    d["date"] = pd.to_datetime(d.date, errors="raise").dt.strftime("%Y-%m-%d")
    if d.date.isna().any() or d.date.duplicated().any():
        raise ValueError(f"{path} 日期有空值或重复；必须一日一行")
    d["tide_height_m"] = pd.to_numeric(d.tide_height_m, errors="raise")
    if np.isinf(d.tide_height_m).any() or (not allow_missing and d.tide_height_m.isna().any()):
        raise ValueError(f"{path} 含未匹配潮位；先明确处理缺测日期后再分析")
    for col in ("tide_datum", "tide_source"):
        if col in d and d[col].dropna().nunique() > 1:
            raise ValueError(f"{path} 混有多个 {col}，请分站点/基准处理")
    return d


def datum_label(d):
    if "tide_datum" in d and d.tide_datum.notna().any():
        return str(d.tide_datum.dropna().iloc[0])
    return "输入潮位表基准（未注明）"


def daily_tides(path, max_span="30min"):
    d = pd.read_csv(path)
    required = {"overpass_utc", "tide_height_m"}
    if not required.issubset(d):
        raise ValueError(f"场景潮位表缺少 {sorted(required - set(d))}")
    d["overpass_utc"] = pd.to_datetime(d.overpass_utc, utc=True, errors="raise")
    if d.empty or d.overpass_utc.isna().any():
        raise ValueError("场景潮位表为空或含空过境时间")
    d["tide_height_m"] = pd.to_numeric(d.tide_height_m, errors="raise")
    if np.isinf(d.tide_height_m).any():
        raise ValueError("场景潮位含无穷值")
    d["date"] = d.overpass_utc.dt.strftime("%Y-%m-%d")
    span = pd.Timedelta(max_span)
    if pd.isna(span) or span <= pd.Timedelta(0):
        raise ValueError("同日过境时间跨度上限必须大于零")
    rows = []
    meta = [c for c in ("tide_source", "tide_model", "tide_datum", "tide_station_file", "tide_station_name",
                        "tide_input_unit", "tide_input_timezone", "tide_datum_offset_m") if c in d]
    for c in meta:
        if d[c].dropna().nunique() > 1:
            raise ValueError(f"场景潮位表混有多个 {c}，不能合为一个站点")
    for date, g in d.groupby("date", sort=True):
        if g.overpass_utc.max() - g.overpass_utc.min() > span:
            raise ValueError(f"{date} 同日过境跨度超过 {max_span}，不宜共用一个日潮位；需按轨道/时刻拆分分析")
        row = {"date": date, "tide_height_m": g.tide_height_m.mean() if g.tide_height_m.notna().all() else np.nan,
               "tide_n_scenes": len(g), "tide_overpass_span_seconds": (g.overpass_utc.max() - g.overpass_utc.min()).total_seconds()}
        row.update({c: g[c].iloc[0] for c in meta})
        rows.append(row)
    return pd.DataFrame(rows)
