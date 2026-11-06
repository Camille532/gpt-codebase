#!/usr/bin/env python3
"""Supplement existing mangrove outputs; never modifies pipeline inputs or assets.

Dependencies: numpy, pandas, scipy. See README.md for estimands and denominators.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform

import numpy as np
import pandas as pd
from scipy import stats

INDICES = ["NDVI", "MNDWI", "EVI", "NDRE", "NIR_B8", "SWIR_B11"]
LABELS = {"core": "核心区", "transition": "原过渡带", "seaward": "向海林缘",
          "landward": "向陆林缘", "mixed": "混合林缘", "whole": "全林",
          "focal": "向海林缘＋过渡带", "rest": "其余森林",
          "transition_seaward": "靠潮汐区一侧的过渡带", "transition_other": "其余过渡带"}
PAIRS = [
    ("rf_density_fixed", "biomass_fixed_mask", "RF_trainLow_imgLow", "RF_trainLow_imgHigh"),
    ("rf_density_retrained", "biomass_fixed_mask", "RF_trainLow_imgLow", "RF_trainHigh_imgHigh"),
    ("linear_ndvi", "biomass_fixed_mask", "LIN_NDVI_imgLow", "LIN_NDVI_imgHigh"),
    ("linear_evi", "biomass_fixed_mask", "LIN_EVI_imgLow", "LIN_EVI_imgHigh"),
    ("linear_nirv", "biomass_fixed_mask", "LIN_NIRv_imgLow", "LIN_NIRv_imgHigh"),
    ("R1", "mapped_area_ha", "R1_NDVI_MNDWImask_imgLow", "R1_NDVI_MNDWImask_imgHigh"),
    ("R2", "mapped_area_ha", "R2_NDVIonly_imgLow", "R2_NDVIonly_imgHigh"),
    ("RF_fixed", "mapped_area_ha", "RFC_trainLow_imgLow", "RFC_trainLow_imgHigh"),
    ("RF_retrained", "mapped_area_ha", "RFC_trainLow_imgLow", "RFC_trainHigh_imgHigh"),
    ("R1", "biomass_in_mapped_area", "R1_NDVI_MNDWImask_imgLow", "R1_NDVI_MNDWImask_imgHigh"),
    ("R2", "biomass_in_mapped_area", "R2_NDVIonly_imgLow", "R2_NDVIonly_imgHigh"),
    ("RF_fixed", "biomass_in_mapped_area", "RFC_trainLow_imgLow", "RFC_trainLow_imgHigh"),
    ("RF_retrained_classifier", "biomass_in_mapped_area", "RFC_trainLow_imgLow", "RFC_trainHigh_imgHigh"),
]


def save(df, path):
    df = df.copy()
    if "part" in df:
        df["part_cn"] = df.part.map(LABELS).fillna(df.part)
    if "zone" in df:
        df["zone_cn"] = df.zone.map(LABELS).fillna(df.zone)
    if "reference" in df:
        df["reference_cn"] = df.reference.map(LABELS).fillna(df.reference)
    df.to_csv(path, index=False, encoding="utf-8-sig")


def read_csv(path, required=()):
    d = pd.read_csv(path)
    missing = set(required) - set(d.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    return d


def date_strings(x):
    return pd.to_datetime(x, errors="raise").dt.strftime("%Y-%m-%d")


def truth(x):
    return x.astype(str).str.strip().str.lower().isin(["true", "1", "1.0"])


def fdr(p):
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    keep = np.flatnonzero(np.isfinite(p))
    order = keep[np.argsort(p[keep])]
    if len(order):
        vals = p[order] * len(order) / np.arange(1, len(order) + 1)
        out[order] = np.minimum.accumulate(vals[::-1])[::-1].clip(0, 1)
    return out


def load_tides(site, override=None):
    path = Path(override) if override else site / "results/1_waterline/date_front_wetness.csv"
    q = read_csv(path, ["date"])
    q["date"] = date_strings(q.date)
    if q.date.duplicated().any():
        raise ValueError("Tide table contains duplicate dates")
    tide_col = "tide_height_m" if "tide_height_m" in q else "tide"
    q["tide"] = pd.to_numeric(q[tide_col], errors="raise")
    if "qualified_strict" in q:
        q["strict"] = truth(q.qualified_strict)
    elif "analysis_date" in q:
        # This archived s11 table was created with the strict-date list.
        q["strict"] = truth(q.analysis_date)
    else:
        q["strict"] = False
    q["month"] = pd.to_datetime(q.date).dt.month
    return q[["date", "tide", "strict", "month"]], path


def cluster_meat(scores, codes):
    sums = np.zeros((int(codes.max()) + 1, scores.shape[1]))
    np.add.at(sums, codes, scores)
    return sums.T @ sums


def two_way_cov(X, residual, dates, blocks, k):
    n = len(X)
    codes = [pd.factorize(dates)[0], pd.factorize(blocks)[0],
             pd.factorize(dates.astype(str) + "|" + blocks.astype(str))[0]]
    sizes = [int(c.max()) + 1 for c in codes]
    if min(sizes) < 2 or n <= k:
        raise ValueError("Insufficient clusters or residual degrees of freedom")
    scores = X * residual[:, None]
    meats = [g / (g - 1) * (n - 1) / (n - k) * cluster_meat(scores, c)
             for c, g in zip(codes, sizes)]
    bread = np.linalg.inv(X.T @ X)
    return bread @ (meats[0] + meats[1] - meats[2]) @ bread


def fit_panel(d, group_names):
    """Exact within-point + FWL date/month controls used in original s13.

    All outcomes must have identical finite rows. Joint projection saves time;
    no observations are pooled or averaged before fitting.
    """
    d = d.copy().reset_index(drop=True)
    kx = len(group_names) - 1
    X = np.column_stack([d.tide.to_numpy() * d.group.eq(g).to_numpy()
                         for g in group_names if g != 0])
    controls = [d.group.eq(g).to_numpy() * d.month.eq(m).to_numpy()
                for g in group_names if g != 0 for m in range(2, 13)]
    Z = np.column_stack(controls + [pd.get_dummies(d.date).to_numpy(float)])
    Y = d[INDICES].to_numpy(float)
    W = np.column_stack([Y, X, Z])
    W -= pd.DataFrame(W).groupby(d.pt_id.to_numpy()).transform("mean").to_numpy()
    A, Z = W[:, :len(INDICES) + kx], W[:, len(INDICES) + kx:]
    coef, _, rank_z, _ = np.linalg.lstsq(Z, A, rcond=None)
    A = A - Z @ coef
    Y, X = A[:, :len(INDICES)], A[:, len(INDICES):]
    if np.linalg.matrix_rank(X) != kx:
        raise ValueError("Tide × zone contrasts are not identified in this subset")
    beta = np.linalg.solve(X.T @ X, X.T @ Y)
    residuals = Y - X @ beta
    k = kx + int(rank_z) + d.pt_id.nunique()
    ordered = [g for g in group_names if g != 0]
    return beta, [two_way_cov(X, residuals[:, j], d.date, d.block_id, k)
                  for j in range(len(INDICES))], ordered


def contrast_rows(d, group_names, model):
    beta, covs, ordered = fit_panel(d, group_names)
    comparisons = [(2, 0), (1, 0), (2, 1)]
    if 5 in group_names:
        comparisons += [(5, 0), (1, 5)]
    dof = min(d.date.nunique(), d.block_id.nunique()) - 1
    rows = []
    for a, b in comparisons:
        v = np.array([float(g == a) - float(g == b) for g in ordered])
        for j, index in enumerate(INDICES):
            estimate = float(v @ beta[:, j])
            variance = float(v @ covs[j] @ v)
            # Two-way finite-sample covariance can be non-PSD. Do not turn a
            # negative contrast variance into a zero SE or an artificial p=0.
            se = np.sqrt(variance) if variance >= 0 else np.nan
            p = 2 * stats.t.sf(abs(estimate / se), dof) if se > 0 else np.nan
            width = stats.t.ppf(.975, dof) * se
            rows.append(dict(model=model, index=index,
                contrast=f"{group_names[a]}_minus_{group_names[b]}",
                zone=group_names[a], reference=group_names[b],
                estimate=estimate, se=se, ci95_low=estimate-width, ci95_high=estimate+width,
                error_minus=width, error_plus=width, p_t_cluster=p,
                variance_status="ok" if variance >= 0 else "negative_two_way_variance",
                n_obs=len(d), n_dates=d.date.nunique(), n_points=d.pt_id.nunique(),
                n_blocks=d.block_id.nunique(), unit="index_or_reflectance_units_per_m"))
    out = pd.DataFrame(rows)
    out["q_fdr_six_indices"] = out.groupby("contrast")["p_t_cluster"].transform(
        lambda x: fdr(x.to_numpy()))
    return out


def spectral_analysis(site, q, out):
    path = site / "inputs/cg_edge_spectra_v4_fullyear.csv"
    d = read_csv(path, ["pt_id", "date", "block_id", "edge_class", "d_tidal", "d_nontidal"])
    d["date"] = date_strings(d.date)
    if d.duplicated(["pt_id", "date"]).any():
        raise ValueError("Duplicate point/date spectral records")
    if d.groupby("pt_id").edge_class.nunique().max() > 1:
        raise ValueError("Point zones change across dates")
    d["NIR_B8"], d["SWIR_B11"] = d.B8 / 1e4, d.B11 / 1e4
    d = d[d.edge_class.isin([0, 1, 2, 3])].merge(q, on="date", validate="many_to_one")
    n_before = len(d)
    d = d[np.isfinite(d[INDICES + ["tide", "d_tidal", "d_nontidal"]]).all(axis=1)].copy()
    if d.empty:
        raise ValueError("No finite spectral records match tide dates")
    tables, profiles = [], []
    for refined in [False, True]:
        names = {0: "core", 1: "transition", 2: "seaward", 3: "landward"}
        z = d.copy()
        z["group"] = z.edge_class.astype(int)
        if refined:
            # Geometry-defined proxy, not proof of tidal access or adjacency.
            z.loc[z.edge_class.eq(1) & z.d_tidal.ge(z.d_nontidal), "group"] = 5
            names[1], names[5] = "transition_seaward", "transition_other"
        for subset in ["all", "strict"]:
            use = z if subset == "all" else z[z.strict]
            if use.date.nunique() < 15:
                print(f"Skip spectral {refined=}/{subset}: fewer than 15 dates", flush=True)
                continue
            model = ("refined_" if refined else "original_") + subset
            print(f"Spectral {model}: {len(use)} records, {use.date.nunique()} dates", flush=True)
            tables.append(contrast_rows(use, names, model))
            for g, label in names.items():
                sub = use[use.group.eq(g)]
                profiles.append(dict(model=model, zone=label, n_points=sub.pt_id.nunique(),
                                     n_obs=len(sub), n_dates=sub.date.nunique()))
    result = pd.concat(tables, ignore_index=True)
    save(result, out / "01_spectral_contrasts.csv")
    save(pd.DataFrame(profiles), out / "01_zone_samples.csv")
    return dict(spectral_dates=int(d.date.nunique()), spectral_records=int(len(d)),
                spectral_joint_complete_rows_dropped=int(n_before-len(d)),
                strict_spectral_dates=int(d.loc[d.strict, "date"].nunique()),
                refined_transition_points=int(d.loc[d.edge_class.eq(1) & d.d_tidal.lt(d.d_nontidal), "pt_id"].nunique()))


def partitions(refined=False):
    parts = {"core": [0], "transition_seaward" if refined else "transition": [1],
             "seaward": [2], "landward": [3], "mixed": [4]}
    if refined:
        parts["transition_other"] = [5]
    parts["whole"] = sorted(g for v in parts.values() for g in v)
    parts["focal"] = [1, 2]
    parts["rest"] = [g for g in parts["whole"] if g not in [1, 2]]
    return parts


def summary_vector(s, metric, scenario, groups):
    z = s[(s.metric == metric) & (s.scenario == scenario) & s.grp.isin(groups)]
    if z.grp.duplicated().any():
        raise ValueError(f"Duplicate summary groups: {metric}/{scenario}")
    if set(z.grp) != set(groups):
        raise ValueError(f"Incomplete groups for {metric}/{scenario}; missing {set(groups)-set(z.grp)}")
    return z.set_index("grp").value.reindex(groups).astype(float)


def composite_analysis(site, out, refined=False):
    suffix = "_local_whole" if refined else ""
    s = read_csv(site / f"inputs/cg_s15_summary{suffix}.csv", ["site", "metric", "scenario", "grp", "value", "ref_name"])
    s["grp"] = pd.to_numeric(s.grp, errors="raise").astype(int)
    s = s[s.grp.ge(0)].copy()
    s["value"] = pd.to_numeric(s.value, errors="raise")
    if s.site.nunique() != 1 or s.ref_name.nunique() != 1:
        raise ValueError("Run one site and one biomass reference at a time")
    is_agb = s.ref_name.iloc[0] == "agb_Mg_per_ha"
    parts = partitions(refined)
    groups = parts["whole"]
    refs = summary_vector(s, "mapped_area_ha", "V4_REFERENCE", groups)
    profile = [dict(part=k, area_ha=float(refs.loc[v].sum()),
                    share_whole_area_pct=100 * refs.loc[v].sum()/refs.sum()) for k, v in parts.items()]
    save(pd.DataFrame(profile), out / "00_zone_area_profile.csv")
    rows, audits = [], []
    for method, metric, lo, hi in PAIRS:
        A, B = summary_vector(s, metric, lo, groups), summary_vector(s, metric, hi, groups)
        variants = [(metric, "ha" if metric == "mapped_area_ha" else ("Mg_AGB" if is_agb else "reference_units_ha"), 1.)]
        if metric == "biomass_in_mapped_area" and is_agb:
            variants.append(("aboveground_carbon_in_mapped_area", "Mg_C", .47))
        for output_metric, unit, factor in variants:
            a, b = A * factor, B * factor
            total_low, total_delta = float(a.sum()), float((b-a).sum())
            sum_abs = float((b-a).abs().sum())
            for part, ids in parts.items():
                low, high = float(a.loc[ids].sum()), float(b.loc[ids].sum())
                delta = high-low
                rows.append(dict(method=method, metric=output_metric, unit=unit, part=part,
                    low=low, high=high, delta=delta,
                    relative_change_pct=100*delta/low if low else np.nan,
                    contribution_to_whole_change_pp=100*delta/total_low if total_low else np.nan,
                    whole_low=total_low, low_scenario=lo, high_scenario=hi))
            focal_delta = float((b-a).loc[parts["focal"]].sum())
            rest_delta = float((b-a).loc[parts["rest"]].sum())
            if not np.isclose(focal_delta + rest_delta, total_delta, rtol=1e-9, atol=1e-6):
                raise AssertionError("Local/whole additivity failed")
            audits.append(dict(method=method, metric=output_metric, unit=unit,
                whole_delta=total_delta, focal_delta=focal_delta, rest_delta=rest_delta,
                sum_absolute_component_deltas=sum_abs,
                cancellation_fraction=1-abs(total_delta)/sum_abs if sum_abs else 0.,
                additivity_error=focal_delta+rest_delta-total_delta))
    save(pd.DataFrame(rows), out / "02_composite_local_whole.csv")
    save(pd.DataFrame(audits), out / "03_contribution_checks.csv")
    return dict(summary_reference=s.ref_name.iloc[0], carbon_conversion_used=.47 if is_agb else None,
                composite_reference_area_ha=float(refs.sum()))


def make_daily(m, q, refined=False, min_cover=.8, core_r1=.8):
    m = m.copy()
    m["grp"] = pd.to_numeric(m.grp, errors="raise").astype(int)
    for col in ["n_valid", "n_R1", "n_R2", "n_RF"]:
        m[col] = pd.to_numeric(m[col], errors="raise")
    parts = partitions(refined)
    groups = parts["whole"]
    ref = m[m.date.eq("V4_TOTAL") & m.grp.isin(groups)]
    if set(ref.grp) != set(groups) or ref.grp.duplicated().any():
        raise ValueError("Missing or duplicate V4_TOTAL groups")
    ref = ref.set_index("grp").n_valid.reindex(groups)
    if ref.le(0).any():
        raise ValueError("Empty reference zone")
    base = {part: float(ref.loc[ids].sum())*.01 for part, ids in parts.items()}
    x = m[~m.date.eq("V4_TOTAL") & m.grp.isin(groups)].copy()
    x["date"] = date_strings(x.date)
    if x.duplicated(["date", "grp"]).any():
        raise ValueError("Duplicate date/group map counts")
    rows = []
    for date, g in x.groupby("date"):
        g = g.set_index("grp").reindex(groups)
        valid = g.n_valid
        if valid.gt(ref * (1+1e-6)).any() or valid.lt(0).any():
            raise ValueError(f"Invalid coverage counts on {date}")
        r = {"date": date, "cover": float(valid.sum()/ref.sum())}
        for part, ids in parts.items():
            r[f"cover_{part}"] = float(valid.loc[ids].sum()/ref.loc[ids].sum())
        for meth in ["R1", "R2", "RF"]:
            counts = g["n_"+meth]
            if counts.lt(0).any() or counts.gt(valid*(1+1e-6)).any():
                raise ValueError(f"Classified pixels exceed valid pixels: {date}/{meth}")
            by_group = counts/valid.where(valid.gt(0)) * ref*.01
            for part, ids in parts.items():
                r[f"A_{meth}_{part}_ha"] = float(by_group.loc[ids].sum(min_count=len(ids)))
        rows.append(r)
    a = pd.DataFrame(rows).merge(q[["date", "tide"]], on="date", validate="one_to_one")
    a["year"], a["month"] = pd.to_datetime(a.date).dt.year, pd.to_datetime(a.date).dt.month
    ok = a.cover.ge(min_cover) & a.A_R1_core_ha.div(base["core"]).ge(core_r1) & np.isfinite(a.tide)
    a["main_qc_pass"] = ok
    a["local_coverage_pass"] = a.cover_seaward.ge(min_cover) & a[
        "cover_transition_seaward" if refined else "cover_transition"].ge(min_cover)
    return a, base


def daily_effects(a, base, model):
    if len(a) < 20:
        return pd.DataFrame()
    X = np.column_stack([np.ones(len(a)), a.tide.to_numpy(float),
        pd.get_dummies(a.year, drop_first=True).to_numpy(float),
        pd.get_dummies(a.month, drop_first=True).to_numpy(float)])
    rows = []
    for meth in ["R1", "R2", "RF"]:
        for part, area in base.items():
            y = a[f"A_{meth}_{part}_ha"].to_numpy(float)/area*100
            ok = np.isfinite(y) & np.isfinite(X).all(axis=1)
            xx, yy = X[ok], y[ok]
            rank = np.linalg.matrix_rank(xx)
            if rank != xx.shape[1] or len(yy) <= rank:
                raise ValueError(f"Unidentified daily-area model: {model}/{meth}/{part}")
            b = np.linalg.lstsq(xx, yy, rcond=None)[0]
            e = yy - xx@b
            bread = np.linalg.inv(xx.T@xx)
            V = bread @ ((xx.T*e**2)@xx) @ bread * len(yy)/(len(yy)-rank)
            se = float(np.sqrt(max(V[1,1], 0)))
            rows.append(dict(model=model, method=meth, part=part, n_dates=len(yy),
                slope_part_pp_per_m=b[1], se_hc1=se,
                ci95_low=b[1]-1.96*se, ci95_high=b[1]+1.96*se,
                p_normal_hc1=2*stats.norm.sf(abs(b[1]/se)) if se else np.nan,
                slope_ha_per_m=b[1]*area/100,
                contribution_whole_pp_per_m=b[1]*area/base["whole"],
                reference_area_ha=area))
    return pd.DataFrame(rows)


def paired_changes(a, base, early_end=2020, late_start=2023):
    early, late = a[a.year.le(early_end)], a[a.year.ge(late_start)]
    rows = []
    for e in early.itertuples(index=False):
        for l in late.itertuples(index=False):
            month_gap = abs(e.month-l.month)
            if min(month_gap, 12-month_gap) > 1:
                continue
            for meth in ["R1", "R2", "RF"]:
                for part, area in base.items():
                    delta = getattr(l, f"A_{meth}_{part}_ha")-getattr(e, f"A_{meth}_{part}_ha")
                    if not np.isfinite(delta):
                        continue
                    rows.append(dict(early_date=e.date, late_date=l.date, method=meth, part=part,
                        dtide_m=l.tide-e.tide, delta_ha=delta,
                        delta_part_reference_pp=100*delta/area,
                        delta_whole_reference_pp=100*delta/base["whole"]))
    P = pd.DataFrame(rows)
    if P.empty:
        return P, pd.DataFrame()
    labels = ["late_lower_gt1m", "late_lower_0.3to1m", "similar_within0.3m",
              "late_higher_0.3to1m", "late_higher_gt1m"]
    P["tide_bin"] = pd.cut(P.dtide_m, [-np.inf,-1.,-.3,.3,1.,np.inf], labels=labels)
    S = P.groupby(["method", "part", "tide_bin"], observed=True).agg(
        n_pairs=("delta_ha", "size"), n_unique_early=("early_date", "nunique"),
        n_unique_late=("late_date", "nunique"), mean_delta_ha=("delta_ha", "mean"),
        mean_part_reference_pp=("delta_part_reference_pp", "mean"),
        mean_whole_reference_pp=("delta_whole_reference_pp", "mean"),
        p10_part_reference_pp=("delta_part_reference_pp", lambda x: x.quantile(.1)),
        p90_part_reference_pp=("delta_part_reference_pp", lambda x: x.quantile(.9))).reset_index()
    return P, S


def daily_analysis(site, q, out, refined=False):
    suffix = "_local_whole" if refined else ""
    maps = read_csv(site / f"inputs/cg_s16_daily_maps{suffix}.csv", ["date", "grp", "n_valid", "n_R1", "n_R2", "n_RF"])
    a, base = make_daily(maps, q, refined)
    save(a, out / "04_daily_area_and_qc.csv")
    main = a[a.main_qc_pass].copy()
    stricter = main[main.local_coverage_pass].copy()
    eff = [daily_effects(main, base, "original_qc")]
    if len(stricter) >= 20:
        eff.append(daily_effects(stricter, base, "focal_coverage_ge80pct"))
    save(pd.concat(eff, ignore_index=True), out / "05_daily_area_effects.csv")
    P, S = paired_changes(main, base)
    save(P, out / "06_pair_details.csv")
    save(S, out / "06_pair_summary.csv")
    return dict(map_export_dates=int(maps.loc[~maps.date.eq("V4_TOTAL"), "date"].nunique()),
                daily_area_dates=int(len(main)), focal_coverage_qc_dates=int(len(stricter)),
                n_unique_date_pairs=int(P[["early_date", "late_date"]].drop_duplicates().shape[0]) if len(P) else 0,
                daily_reference_area_ha=base["whole"])


def factorial_analysis(site, out, refined):
    path = site / "inputs/cg_s17_factorial_common.csv"
    if not path.exists():
        return {"factorial_available": False}
    s = read_csv(path, ["method", "scenario", "grp", "value", "ref_name"])
    s["grp"] = pd.to_numeric(s.grp).astype(int)
    s["value"] = pd.to_numeric(s.value)
    groups = partitions(refined)
    if (5 in s.grp.unique()) != refined:
        print("Skip factorial export: partition differs from current run", flush=True)
        return {"factorial_available": False, "factorial_skip_reason": "different_partition"}
    rows = []
    for method in ["R1", "R2", "RF"]:
        z = s[s.method.eq(method)]
        if z.ref_name.nunique() != 1 or z.ref_name.iloc[0] != "agb_Mg_per_ha":
            raise ValueError("Factorial carbon decomposition requires an AGB reference")
        wide = z.pivot(index="grp", columns="scenario", values="value")
        for part, ids in groups.items():
            t = wide.loc[ids].sum(min_count=len(ids))
            a, b, c, d = [float(t[k])*.47 for k in ["Mlow_Dlow", "Mhigh_Dlow", "Mlow_Dhigh", "Mhigh_Dhigh"]]
            mask = .5*((b-a)+(d-c))
            density = .5*((c-a)+(d-b))
            if not np.isclose(mask+density, d-a, atol=1e-6):
                raise AssertionError("Factorial additivity failure")
            rows.append(dict(method=method, part=part, unit="Mg_C", low=a, high=d,
                delta=d-a, mask_component=mask, density_component=density,
                interpretation="symmetric_model_output_decomposition_not_causal_attribution"))
    save(pd.DataFrame(rows), out / "07_carbon_factorial_decomposition.csv")
    coverage = s[s.method.eq("COVERAGE")].pivot(index="grp", columns="scenario", values="value")
    coverage["common_fraction"] = coverage.common_ha/coverage.reference_ha
    save(coverage.reset_index(), out / "07_common_footprint_coverage.csv")
    return {"factorial_available": True}


def write_readout(out, meta):
    lines = ["# 本次补充分析说明", "",
        "这些结果评估观测与模型输出的敏感性，不是独立精度验证，也不是潮位的因果效应。", "",
        "## 日期与版本", "", "```json", json.dumps(meta, ensure_ascii=False, indent=2), "```", "",
        "## 先看哪些表", "",
        "1. 01_spectral_contrasts.csv：original_all为原分区；refined_all拆分过渡带。seaward_minus_transition是直接差异检验。",
        "2. 02_composite_local_whole.csv：relative_change_pct用各分区自己的低潮结果作分母；contribution_to_whole_change_pp统一用全林低潮结果作分母，可相加。",
        "3. 03_contribution_checks.csv：focal_delta＋rest_delta＝whole_delta；cancellation_fraction越接近1，分区正负变化抵消越多。",
        "4. 05_daily_area_effects.csv：过渡带、向海林缘、全林的逐日面积斜率；original_qc复用旧质控，另提供局部覆盖率敏感性分析。",
        "5. 06_pair_summary.csv：前后期月份接近的配对，百分比以固定分区参考面积为分母；配对重复使用日期，不视作独立样本。", "",
        "## 解释限制", "",
        "- 原过渡带是所有边界内侧30–60 m，不是已划定的向海过渡带或潮间带。",
        "- 新过渡带用d_tidal < d_nontidal定义，是几何代理，不能证明淹没或水流连通，也不是严格沿同一横断面的距离梯度。",
        "- 光谱模型有点、日期固定效应及分区×月份控制，只识别相对核心区的斜率差。",
        "- 光谱补充采用日期×空间区块双向聚类协方差、较小簇数减1的t检验与t区间；未重跑旧WCR。不能把本表p值标成WCR结果。",
        "- FDR在每个模型、每个对比的6个指标内调整；多个模型是补充/敏感性分析，不是全部检验共同校正。",
        "- 逐日面积回归沿用HC1及年份、月份控制，未处理全部时间自相关，年份控制也不等于已验证真实森林变化。",
        "- 生物量参考是产品而非本地独立实测；碳量仅在AGB参考下乘0.47，为对应地上碳。",
        "- 原高低潮合成汇总未强制共同有效像元范围，不能解释成纯潮位效应；可选GEE脚本另导出共同范围的碳量分解。",
        "- 合成差异不等于精度误差，前后期面积差异不等于已证实的虚假损失。",
        "- 原GEE高低潮组按月平衡，未完全平衡年份、天气、观测数量。",
        "- 本次读取已有潮位和导出结果，没有重新运行未移时FES；不将旧FES结果改称0小时FES结果。",
        "- g08按pixelArea汇总，g09沿用像元数×0.01ha，两类参考面积略有差异，不交叉混用。", ""]
    (out / "READ_RESULTS.md").write_text("\n".join(lines), encoding="utf-8")


def date_audit(site, out):
    path = site / "results/1_waterline/date_front_wetness.csv"
    if not path.exists():
        return {}
    w = read_csv(path, ["date", "tide", "wet_frac_bin1", "n_valid_bin1"])
    valid = w.wet_frac_bin1.notna() & w.n_valid_bin1.gt(0) & np.isfinite(w.tide)
    strict = truth(w.analysis_date) if "analysis_date" in w else pd.Series(False, index=w.index)
    meta = dict(waterline_tide_table_dates=int(len(w)),
        front_valid_optical_dates=int(valid.sum()),
        front_majority_wet_dates=int((valid & w.wet_frac_bin1.gt(.5)).sum()),
        front_strict_valid_dates=int((valid & strict).sum()))
    save(pd.DataFrame([dict(quantity=k, count=v) for k,v in meta.items()]), out/"00_date_audit.csv")
    return meta


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--site-dir", type=Path, required=True, help="Folder containing inputs/ and results/")
    ap.add_argument("--out-dir", type=Path)
    ap.add_argument("--dates", type=Path, help="Optional original dates CSV; otherwise use archived waterline tide table")
    ap.add_argument("--skip-spectra", action="store_true")
    ap.add_argument("--refined-exports", action="store_true", help="Read new GEE summary/maps with group5 = other transition")
    args = ap.parse_args()
    site = args.site_dir.resolve()
    out = (args.out_dir or site / ("results/6_local_whole_refined" if args.refined_exports else "results/5_local_whole")).resolve()
    if out == site or out == site/"inputs" or site/"inputs" in out.parents:
        raise ValueError("Output directory must not overwrite inputs")
    out.mkdir(parents=True, exist_ok=True)
    q, tide_path = load_tides(site, args.dates)
    meta = dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
        site_dir=str(site), partition="refined" if args.refined_exports else "original",
        tide_file=str(tide_path), tide_table_dates=int(q.date.nunique()),
        spectral_inference="two_way_cluster_covariance_with_t_test_no_wild_bootstrap",
        strict_flag_source="qualified_strict" if args.dates else "archived_s11_analysis_date",
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    files = [tide_path] + sorted((site/"inputs").glob("*.csv"))
    meta["input_sha256"] = {str(p.relative_to(site)) if p.is_relative_to(site) else str(p):
        hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    meta.update(date_audit(site, out))
    if not args.skip_spectra:
        meta.update(spectral_analysis(site, q, out))
    print("Summarize composite contrasts and regional contributions", flush=True)
    meta.update(composite_analysis(site, out, args.refined_exports))
    print("Rebuild daily areas and extend paired changes to local zones", flush=True)
    meta.update(daily_analysis(site, q, out, args.refined_exports))
    meta.update(factorial_analysis(site, out, args.refined_exports))
    (out/"run_meta.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
    write_readout(out, meta)
    print(json.dumps({k:v for k,v in meta.items() if k != "input_sha256"}, ensure_ascii=False, indent=2))
    print("DONE:", out)


if __name__ == "__main__":
    main()
