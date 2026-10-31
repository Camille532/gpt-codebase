#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Write separate GEE script copies with observed date groups and a rerun audit."""
from pathlib import Path
import argparse
import json
import re
from tide_data import read_date_table, datum_label
from make_tide_composite_dates import select_groups


def array_pattern(name):
    return re.compile(r"\bvar\s+" + re.escape(name) + r"\s*=\s*\[.*?\]\s*;", re.S)


def read_array(text, name):
    hits = list(array_pattern(name).finditer(text))
    if len(hits) != 1:
        raise ValueError(f"应有且仅有一个var {name}数组，实际{len(hits)}个")
    return re.findall(r"['\"](\d{4}-\d{2}-\d{2})['\"]", hits[0].group())


def replace_array(text, name, values):
    before = read_array(text, name)
    rows = [", ".join(repr(d) for d in values[i:i + 8]) for i in range(0, len(values), 8)]
    block = f"var {name} = [\n  " + ",\n  ".join(rows) + "\n];"
    text = array_pattern(name).sub(lambda _: block, text)
    return text, {"old_count": len(before), "new_count": len(values),
                  "added": sorted(set(values) - set(before)), "removed": sorted(set(before) - set(values)),
                  "changed": set(before) != set(values)}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dates", type=Path, required=True, help="更新后的合格分析日期表")
    ap.add_argument("--waterline-dates", type=Path, help="水边线用日期表；不传则保留g06原DATES")
    ap.add_argument("--gee-dir", type=Path, required=True, help="包含cg_g06至g09原脚本的目录")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()
    if args.gee_dir.resolve() == args.out_dir.resolve():
        raise ValueError("输出目录须与原GEE目录分开")
    d = read_date_table(args.dates)
    if "tide_source" not in d or not d.tide_source.eq("observed_station").all():
        raise ValueError("日期表必须由实测潮位更新生成")
    dates = sorted(d.date.tolist())
    hi, lo = select_groups(d)
    groups = {"cg_g06_waterline.js": {}, "cg_g07_spectra.js": {"SELECTED_DATES": dates},
              "cg_g08_tide_composite_experiment.js": {"HIGH_DATES": hi, "LOW_DATES": lo},
              "cg_g09_change_detection_sim.js": {"DATES": dates, "LOW_DATES": lo}}
    if args.waterline_dates:
        groups["cg_g06_waterline.js"]["DATES"] = sorted(read_date_table(args.waterline_dates).date.tolist())
    audit, staged = {}, []
    for name, arrays in groups.items():
        path = args.gee_dir / name
        text = path.read_text(encoding="utf-8")
        audit[name] = {}
        for variable, values in arrays.items():
            text, audit[name][variable] = replace_array(text, variable, values)
        text = re.sub(r"^// 高潮组均值.*$", "// 潮位均值应以本次实测分组输出为准。", text, flags=re.M)
        text = re.sub(r"var TIDE_DATE_SOURCE = ['\"]legacy_fes['\"];", "var TIDE_DATE_SOURCE = 'observed_station';", text)
        text = "// 实测日期更新；潮位基准：" + datum_label(d) + "；输入：" + args.dates.name + "\n" + text
        staged.append((name, text))
    g6 = audit["cg_g06_waterline.js"].get("DATES", {"added": []})
    g7 = audit["cg_g07_spectra.js"]["SELECTED_DATES"]
    g8 = audit["cg_g08_tide_composite_experiment.js"]
    g9 = audit["cg_g09_change_detection_sim.js"]
    plan = {"g06_export_new_dates_only": g6["added"], "g07_export_new_dates_only": g7["added"],
            "g08_rerun": g8["HIGH_DATES"]["changed"] or g8["LOW_DATES"]["changed"],
            "g09_R1_R2_export_new_dates_only": g9["DATES"]["added"],
            "g09_RF_retrain_and_export_all_selected_dates": g9["LOW_DATES"]["changed"],
            "g09_RF_export_new_dates_only_if_training_unchanged": g9["DATES"]["added"],
            "local_recompute": ["s11", "s13 (或已使用的s08/s09)", "s16", "s15 (若g08重跑)"]}
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for name, text in staged:
        (args.out_dir / name).write_text(text, encoding="utf-8")
    (args.out_dir / "rerun_plan.json").write_text(json.dumps({"arrays": audit, "rerun": plan}, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(plan, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
