#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""s15准备：按“同月配对”挑选高潮/低潮日期，用于构建季节构成相同的高潮合成与低潮合成。

规则：在每个日历月内（多年合并）按过境潮位排序，取上三分之一为高潮组、下三分之一为低潮组；
该月合格日期少于3个则跳过。这样两组的月份分布相同，合成影像差异主要来自潮位而非物候。
输出：打印两组日期的潮位与月份统计，并写出GEE可粘贴的 JS 片段。
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from tide_data import read_date_table, datum_label


def select_groups(d):
    d = d.copy()
    d["month"] = pd.to_datetime(d.date).dt.month
    hi, lo = [], []
    for month, g in d.groupby("month"):
        if len(g) < 3:
            continue
        g = g.sort_values(["tide_height_m", "date"], kind="stable")
        k = len(g) // 3
        low, high = g.iloc[:k], g.iloc[-k:]
        # Rounded/equal levels must not manufacture a tide contrast.
        if low.tide_height_m.max() >= high.tide_height_m.min():
            print(f"跳过{month}月：高低潮边界重叠或潮位相同")
            continue
        lo += low.date.tolist()
        hi += high.date.tolist()
    if not hi or not lo:
        raise ValueError("没有可形成潮位差的同月高低潮分组")
    return sorted(hi), sorted(lo)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dates", type=Path, required=True, help="含 date、tide_height_m 的合格日期表")
    ap.add_argument("--out-js", type=Path, required=True)
    ap.add_argument("--site", required=True, help="写入JS注释的站点名")
    args = ap.parse_args()
    d = read_date_table(args.dates)
    d["month"] = pd.to_datetime(d.date).dt.month
    hi, lo = select_groups(d)
    H, L = d[d.date.isin(hi)], d[d.date.isin(lo)]
    print(f"{args.site}: 高潮组 {len(H)} 日，潮位均值 {H.tide_height_m.mean():+.2f} m（{H.tide_height_m.min():+.2f}~{H.tide_height_m.max():+.2f}）；"
          f"低潮组 {len(L)} 日，均值 {L.tide_height_m.mean():+.2f} m（{L.tide_height_m.min():+.2f}~{L.tide_height_m.max():+.2f}）")
    print("两组月份分布：", pd.concat([H.month.value_counts().rename("高"), L.month.value_counts().rename("低")], axis=1)
          .fillna(0).astype(int).sort_index().T.to_string())

    def js(name, dates):
        dates = sorted(dates)
        rows = [",".join(f"'{x}'" for x in dates[i:i + 8]) for i in range(0, len(dates), 8)]
        return f"var {name} = [\n " + ",\n ".join(rows) + "];\n"
    args.out_js.parent.mkdir(parents=True, exist_ok=True)
    args.out_js.write_text(
        f"// {args.site}：同月配对高潮/低潮日期（每月上/下三分之一），由 make_tide_composite_dates.py 生成\n"
        f"// 高潮组均值 {H.tide_height_m.mean():+.2f} m，低潮组均值 {L.tide_height_m.mean():+.2f} m（基准：{datum_label(d)}）\n"
        + js("HIGH_DATES", hi) + js("LOW_DATES", lo), encoding="utf-8")
    print("写出", args.out_js)


if __name__ == "__main__":
    main()
