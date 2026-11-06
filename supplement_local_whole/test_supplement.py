"""Meaningful checks for contrast identification, conservation, and safe JS edits."""
import unittest
import numpy as np
import pandas as pd

from local_whole import INDICES, fit_panel, summary_vector, paired_changes
from make_supplement_gee import transform


class SupplementTests(unittest.TestCase):
    def test_direct_contrast_matches_change_of_reference(self):
        rng = np.random.default_rng(20261106)
        dates = pd.date_range("2019-01-01", periods=36, freq="10D")
        rows = []
        true = [0., .01, .04, -.02]
        day_effects = rng.normal(size=len(dates))
        for group in range(4):
            for point in range(8):
                point_effect = rng.normal()
                for t, date in enumerate(dates):
                    tide = np.sin(t*.71)
                    value = point_effect + day_effects[t] + true[group]*tide + .001*rng.normal()
                    rows.append(dict(pt_id=group*8+point, date=str(date.date()), month=date.month,
                        block_id=group*4+point//2, group=group, tide=tide, **{k:value for k in INDICES}))
        d = pd.DataFrame(rows)
        b, v, order = fit_panel(d, {0:"core", 1:"transition", 2:"seaward", 3:"landward"})
        c = np.array([-1.,1.,0.])
        self.assertAlmostEqual(float(c@b[:,0]), .03, delta=.0005)
        swapped = d.copy()
        swapped["group"] = swapped.group.map({0:1,1:0,2:2,3:3})
        b2, v2, _ = fit_panel(swapped, {0:"transition",1:"core",2:"seaward",3:"landward"})
        self.assertAlmostEqual(float(c@b[:,0]), float(b2[1,0]), places=10)
        self.assertAlmostEqual(float(c@v[0]@c), float(v2[0][1,1]), places=10)

    def test_missing_summary_group_is_not_silently_zero(self):
        s = pd.DataFrame(dict(metric=["m"], scenario=["s"], grp=[0], value=[10.]))
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            summary_vector(s, "m", "s", [0,1])

    def test_pair_denominators_and_season_wrap(self):
        rows = []
        for date, year, month, tide, area in [("2019-12-01",2019,12,0.,100.),
                                             ("2023-01-01",2023,1,1.5,90.)]:
            r = dict(date=date, year=year, month=month, tide=tide)
            r.update({f"A_{m}_{p}_ha": area if p=="whole" else area/2
                      for m in ["R1","R2","RF"] for p in ["whole","seaward"]})
            rows.append(r)
        P, S = paired_changes(pd.DataFrame(rows), {"whole":200.,"seaward":50.})
        row = P[(P.method=="R1") & (P.part=="seaward")].iloc[0]
        self.assertEqual(row.delta_part_reference_pp, -10.)
        self.assertEqual(row.delta_whole_reference_pp, -2.5)
        self.assertEqual(str(row.tide_bin), "late_higher_gt1m")

    def test_generator_rejects_unknown_or_template_sources(self):
        with self.assertRaises(ValueError):
            transform("var ASSET_ROOT='__ASSET_ROOT__';", "g09")
        with self.assertRaises(ValueError):
            transform("var grp = customGrouping();", "g09")
        original = "var grp = edgeAll.rename('grp');\nvar name = PREFIX + 's16_daily_maps';"
        changed = transform(original, "g09")
        self.assertIn("s16_daily_maps_local_whole", changed)
        self.assertIn("gte(cls.select('d_nontidal'))", changed)
        self.assertNotIn("where", original)


if __name__ == "__main__":
    unittest.main(verbosity=2)
