#!/usr/bin/env python3
"""Generate optional GEE scripts beside existing g08/g09; originals untouched.

Only grouping and extra exports change. Dates, features, training points,
classification rules and models remain those in the supplied source scripts.
No GEE assets are created by this Python program. Review and run JS manually.
"""
from pathlib import Path
import argparse
import re
import csv


NEW_GROUP = """// Supplement: group1 = transition nearer tidal zone; group5 = other transition.
// Geometric proxy only; not an intertidal boundary or a hydrological test.
var grp = edgeAll.where(edgeAll.eq(1).and(
  cls.select('d_tidal').gte(cls.select('d_nontidal'))), 5).rename('grp');"""

EXTRA = """
// ----- Optional 2x2 carbon-output decomposition on ONE common valid footprint -----
// The low-trained RF regression is held fixed. All methods use the same density
// images; masks vary between low/high composites. No independent accuracy claim.
var commonS17 = lo.select(FEATS).mask().reduce(ee.Reducer.min())
  .and(hi.select(FEATS).mask().reduce(ee.Reducer.min()))
  .and(preds['RF_trainLow_imgLow'].mask()).and(preds['RF_trainLow_imgHigh'].mask())
  .and(isM).and(zone);
function s17Grouped(img, method, scenario) {
  var r = img.rename('v').addBands(grp).updateMask(zone).reduceRegion({
    reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'grp'}),
    geometry: aoi, crs: PROJ, maxPixels: 1e13, tileScale: 16});
  return ee.List(r.get('groups')).map(function(g) {
    g = ee.Dictionary(g);
    return ee.Feature(null, {method: method, scenario: scenario,
      grp: g.get('grp'), value: g.get('sum'), ref_name: refName});
  });
}
var s17rows = ee.List([]);
[['R1','R1_NDVI_MNDWImask'], ['R2','R2_NDVIonly'], ['RF','RFC_trainLow']].forEach(function(k) {
  var ml = maps[k[1] + '_imgLow'].eq(1).unmask(0);
  var mh = maps[k[1] + '_imgHigh'].eq(1).unmask(0);
  var dl = preds['RF_trainLow_imgLow'], dh = preds['RF_trainLow_imgHigh'];
  [['Mlow_Dlow',ml,dl], ['Mhigh_Dlow',mh,dl],
   ['Mlow_Dhigh',ml,dh], ['Mhigh_Dhigh',mh,dh]].forEach(function(c) {
    s17rows = s17rows.cat(s17Grouped(c[1].multiply(c[2]).multiply(haImg)
      .updateMask(commonS17), k[0], c[0]));
  });
});
s17rows = s17rows.cat(s17Grouped(haImg.updateMask(commonS17), 'COVERAGE', 'common_ha'));
s17rows = s17rows.cat(s17Grouped(haImg.updateMask(isM), 'COVERAGE', 'reference_ha'));
Export.table.toDrive({collection: ee.FeatureCollection(s17rows.flatten()),
  description: PREFIX + 's17_factorial_common', fileFormat: 'CSV',
  selectors: ['method','scenario','grp','value','ref_name']});
"""


def transform(text, kind):
    if re.search(r"__[A-Z_]+__", text):
        raise ValueError("Input is an unrendered template. Use sites/cangio/gee/*.js containing actual dates/assets.")
    pattern = r"var\s+grp\s*=\s*edgeAll\.rename\(['\"]grp['\"]\);"
    text, n = re.subn(pattern, NEW_GROUP, text)
    if n != 1:
        raise ValueError(f"Expected exactly one original grouping statement; found {n}. Source may be a different version.")
    replacements = ["s15_summary", "s15_pixels"] if kind == "g08" else ["s16_daily_maps"]
    for name in replacements:
        old = "'" + name + "'"
        if old not in text:
            raise ValueError(f"Export name not found: {name}")
        text = text.replace(old, "'" + name + "_local_whole'")
    if kind == "g08":
        for token in ["var preds", "var maps", "var FEATS", "var haImg", "var hi", "refName"]:
            if token not in text:
                raise ValueError(f"Unsupported g08 source: missing {token}")
        text += EXTRA
    return "// SUPPLEMENT: original inputs and assets preserved; group1/5 definitions changed.\n" + text


def find_source(folder, needle):
    matches = [p for p in folder.glob("*.js") if needle in p.name and "local_whole" not in p.name]
    if len(matches) != 1:
        raise ValueError(f"Need exactly one {needle} script under {folder}, found {len(matches)}")
    return matches[0]


def dates_in_script(text, name):
    match = re.search(r"var\s+" + re.escape(name) + r"\s*=\s*\[([\s\S]*?)\];", text)
    if not match:
        raise ValueError(f"Cannot read date array {name} from rendered source")
    dates = re.findall(r"['\"](\d{4}-\d{2}-\d{2})['\"]", match.group(1))
    if not dates or len(set(dates)) != len(dates):
        raise ValueError(f"Empty or duplicated dates in {name}")
    return set(dates)


def check_date_versions(site, g08, g09):
    main = dates_in_script(g09, "DATES")
    low = dates_in_script(g09, "LOW_DATES")
    if low != dates_in_script(g08, "LOW_DATES"):
        raise ValueError("g08 and g09 use different low-tide training dates; use the same source version")
    if not (low | dates_in_script(g08, "HIGH_DATES")).issubset(main):
        raise ValueError("Composite dates are outside the g09 main-date set")
    candidates = list((site / "inputs").glob("*s16_daily_maps.csv"))
    if len(candidates) != 1:
        raise ValueError("Need exactly one original *s16_daily_maps.csv to check source-date consistency")
    with candidates[0].open(encoding="utf-8-sig", newline="") as f:
        existing = {r["date"] for r in csv.DictReader(f) if r["date"] != "V4_TOTAL"}
    if main != existing:
        raise ValueError(f"Date-version mismatch: g09 has {len(main)} dates, existing export has {len(existing)}. "
                         "Do not mix observed-tide and FES source scripts/exports.")
    print(f"Date check passed: {len(main)} main dates, {len(low)} low-training dates")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--site-dir", type=Path, required=True)
    ap.add_argument("--g08", type=Path, help="Optional explicit rendered g08 source")
    ap.add_argument("--g09", type=Path, help="Optional explicit rendered g09 source")
    args = ap.parse_args()
    folder = args.site_dir / "gee"
    sources = [("g08", args.g08 or find_source(folder, "g08")),
               ("g09", args.g09 or find_source(folder, "g09"))]
    texts = {k:p.read_text(encoding="utf-8-sig") for k,p in sources}
    check_date_versions(args.site_dir, texts["g08"], texts["g09"])
    # Validate both before writing either output.
    jobs = [(p.with_name(p.stem + "_local_whole.js"), transform(texts[k], k))
            for k, p in sources]
    for path, text in jobs:
        path.write_text(text, encoding="utf-8")
        print("CREATED:", path)
    print("Run both new JS files in GEE and submit Tasks. Download NEW named exports into inputs/.")
    print("Then run: python supplement_local_whole/local_whole.py --site-dir sites/cangio --refined-exports --skip-spectra")


if __name__ == "__main__":
    main()
