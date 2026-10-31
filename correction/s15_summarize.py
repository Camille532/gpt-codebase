#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""s15汇总：高潮合成相对低潮合成，制图面积、生物量与碳储量的偏差（%）。

输入：GEE导出的 hhb_s15_summary.csv、sk_s15_summary.csv（放在 tidedata/s15_composite/）。
输出：tidedata/s15_composite/S15_RESULTS.md 与 s15_bias_table.csv。
碳储量 = 生物量 × 0.47（IPCC默认碳含量），相对偏差与生物量相同。
"""
from pathlib import Path
import argparse
import pandas as pd

CLS = {0: "核心", 1: "过渡", 2: "向海林缘", 3: "向陆林缘", 4: "混合林缘", 255: "林外(评估区)"}
PAIRS = [  # (指标, 低潮情景, 高潮情景, 说明)
    ("biomass_fixed_mask", "RF_trainLow_imgLow", "RF_trainLow_imgHigh", "RF生物量，模型固定，换高潮影像"),
    ("biomass_fixed_mask", "RF_trainLow_imgLow", "RF_trainHigh_imgHigh", "RF生物量，全流程用高潮影像"),
    ("biomass_fixed_mask", "LIN_NDVI_imgLow", "LIN_NDVI_imgHigh", "NDVI线性生物量模型"),
    ("biomass_fixed_mask", "LIN_EVI_imgLow", "LIN_EVI_imgHigh", "EVI线性生物量模型"),
    ("biomass_fixed_mask", "LIN_NIRv_imgLow", "LIN_NIRv_imgHigh", "NIRv线性模型"),
    ("mapped_area_ha", "R1_NDVI_MNDWImask_imgLow", "R1_NDVI_MNDWImask_imgHigh", "制图：NDVI≥0.4且MNDWI<0"),
    ("mapped_area_ha", "R2_NDVIonly_imgLow", "R2_NDVIonly_imgHigh", "制图：仅NDVI≥0.4"),
    ("mapped_area_ha", "RFC_trainLow_imgLow", "RFC_trainLow_imgHigh", "制图：RF分类，模型固定"),
    ("mapped_area_ha", "RFC_trainLow_imgLow", "RFC_trainHigh_imgHigh", "制图：RF分类，全流程高潮"),
    ("biomass_in_mapped_area", "R1_NDVI_MNDWImask_imgLow", "R1_NDVI_MNDWImask_imgHigh", "碳储量：规则①制图×RF生物量"),
    ("biomass_in_mapped_area", "RFC_trainLow_imgLow", "RFC_trainLow_imgHigh", "碳储量：RF制图×RF生物量（模型固定）"),
]


def main():
    root = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path, default=root / "tidedata" / "s15_composite")
    args = ap.parse_args()
    frames = [pd.read_csv(f) for f in sorted(args.dir.glob("*s15_summary*.csv"))]
    if not frames:
        raise SystemExit(f"{args.dir} 中没有 *s15_summary*.csv")
    s = pd.concat(frames, ignore_index=True)
    s = s[pd.to_numeric(s.grp, errors="coerce") >= 0].copy()
    s["grp"] = s.grp.astype(float).astype(int)
    s["value"] = pd.to_numeric(s.value, errors="coerce")
    rows = []
    for site, g in s.groupby("site"):
        ref_name = g.ref_name.dropna().iloc[0]
        v4 = g[(g.scenario == "V4_REFERENCE")].set_index("grp").value
        for metric, a, b, label in PAIRS:
            A = g[(g.metric == metric) & (g.scenario == a)].set_index("grp").value
            B = g[(g.metric == metric) & (g.scenario == b)].set_index("grp").value
            if A.empty or B.empty:
                continue
            mangrove = [k for k in A.index if k != 255]
            # 主口径：v4红树林范围内（0–4类）；林外误判(255)单独报告，避免NDVI阈值把陆生植被计入后掩盖红树林本身的变化
            tot_a, tot_b = A[mangrove].sum(), B[mangrove].sum()
            r = {"site": site, "item": label, "metric": metric, "ref": ref_name,
                 "low_total": tot_a, "high_total": tot_b, "bias_total_pct": 100 * (tot_b - tot_a) / tot_a}
            if 255 in A.index and A[255] > 0:
                r["outside_low"], r["outside_high"] = A[255], B.get(255, 0)
                r["outside_change_pct"] = 100 * (B.get(255, 0) - A[255]) / A[255]
            for k in [0, 1, 2, 3, 4]:
                if k in A.index and A[k] > 0:
                    r[f"bias_{CLS[k]}_pct"] = 100 * (B.get(k, 0) - A[k]) / A[k]
            if metric == "mapped_area_ha":
                r["v4_area_ha"] = v4[[k for k in v4.index if k != 255]].sum()
                r["low_vs_v4_pct"] = 100 * (tot_a - r["v4_area_ha"]) / r["v4_area_ha"]
            rows.append(r)
    t = pd.DataFrame(rows)
    t.to_csv(args.dir / "s15_bias_table.csv", index=False, encoding="utf-8-sig")
    cols = ["site", "item", "bias_total_pct"] + [c for c in t.columns if c.startswith("bias_") and c != "bias_total_pct"] + ["outside_change_pct"]
    md = ["# s15：高潮合成相对低潮合成的偏差（%）", "",
          f"参考变量：{', '.join(sorted(t.ref.unique()))}。碳储量=生物量×0.47，相对偏差与生物量相同。",
          "“全林”= v4红树林范围内（核心+过渡+向海+向陆+混合）；制图结果在林外300 m评估区内的误判面积单独列在 s15_bias_table.csv（outside_*）。", ""]
    hdr = "| " + " | ".join(["站点", "项目", "全林"] + [c.replace("bias_", "").replace("outside_change", "林外误判变化").replace("_pct", "") for c in cols[3:]]) + " |"
    md += [hdr, "|" + "---|" * len(cols)]
    for _, r in t.iterrows():
        md.append("| " + " | ".join([r["site"], r["item"]] + [("" if pd.isna(r[c]) else f"{r[c]:+.1f}") for c in cols[2:]]) + " |")
    (args.dir / "S15_RESULTS.md").write_text("\n".join(md), encoding="utf-8")
    pd.set_option("display.width", 250)
    print(t[cols].round(1).to_string(index=False))


if __name__ == "__main__":
    main()
