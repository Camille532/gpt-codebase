import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
import tempfile
import subprocess
import json
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tide_data import station_table, match_station, utc_times, daily_tides, read_date_table
from refresh_observed_dates import refresh
from make_tide_composite_dates import select_groups
from sync_gee_tide_dates import replace_array
from s11_waterline_elevation import waterline_threshold

class TideMigrationTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.root=Path(self.temp.name)
    def tearDown(self):
        self.temp.cleanup()
    def args(self,path,**kwargs):
        d=dict(station_file=path,station_format='ymdh',station_timezone='UTC',station_unit='mm',
               datum_offset_m=0.,station_sep=',',station_encoding='utf-8-sig',missing_value=[-32767],
               station_time_col='time',station_tide_col='tide_height_m',station_time_format=None,station_sheet=None)
        d.update(kwargs)
        return SimpleNamespace(**d)
    def test_uhslc_no_header_units_and_sentinel_gap(self):
        p=self.root/'raw.csv';p.write_text('2024,1,1,0,1000\n2024,1,1,1,2000\n2024,1,1,2,-32767\n2024,1,1,3,4000\n')
        obs=station_table(self.args(p))
        np.testing.assert_allclose(obs.level,[1,2,4])
        t=['2023-12-31T23:59Z','2024-01-01T00:30Z','2024-01-01T02:00Z','2024-01-01T03:00Z','2024-01-01T03:01Z']
        r=match_station(t,obs,'90min')
        np.testing.assert_allclose(r.tide_height_m,[np.nan,1.5,np.nan,4,np.nan],equal_nan=True)
        self.assertEqual(r.tide_match_status.tolist(),['outside_range','interpolated','gap_too_large','exact','outside_range'])
    def test_timezone(self):
        self.assertEqual(utc_times(['2024-01-01 08:30'],'Asia/Shanghai')[0],pd.Timestamp('2024-01-01T00:30Z'))
        with self.assertRaises(ValueError):utc_times(['2024-01-01 08:30'])
    def test_conflicting_duplicate(self):
        p=self.root/'dup.csv';p.write_text('2024,1,1,0,1000\n2024,1,1,0,2000\n')
        with self.assertRaises(ValueError):station_table(self.args(p))
    def test_exact_in_single_record(self):
        obs=pd.DataFrame({'time':pd.to_datetime(['2024-01-01T00:00Z']),'level':[2.]})
        r=match_station(['2024-01-01T00:00Z','2024-01-01T00:01Z'],obs,'90min')
        self.assertEqual(r.tide_height_m.iloc[0],2.)
        self.assertTrue(pd.isna(r.tide_height_m.iloc[1]))
    def test_daily_missing_is_not_averaged_away(self):
        p=self.root/'sc.csv'
        pd.DataFrame({'overpass_utc':['2024-01-01T03:00Z','2024-01-01T03:01Z'],'tide_height_m':[1,np.nan]}).to_csv(p,index=False)
        self.assertTrue(pd.isna(daily_tides(p).tide_height_m.iloc[0]))
        pd.DataFrame({'overpass_utc':['2024-01-01T03:00Z','2024-01-01T05:00Z'],'tide_height_m':[1,2]}).to_csv(p,index=False)
        with self.assertRaises(ValueError):daily_tides(p)
    def test_refresh_preserves_sample_flags_and_removes_fes(self):
        p=self.root/'old.csv'
        old=pd.DataFrame({'date':['2024-01-01','2024-01-02'],'tide_height_m':[-1,0],'qualified_strict':[True,False],
                          'valid_mean':[.8,.6],'tide_lag_+2.0h':[0,1]})
        old.to_csv(p,index=False)
        new=pd.DataFrame({'date':old.date,'tide_height_m':[3.,2.], 'tide_source':'observed_station'})
        r,a=refresh(p,new)
        self.assertEqual(r.date.tolist(),old.date.tolist())
        self.assertEqual(r.qualified_strict.tolist(),old.qualified_strict.tolist())
        self.assertEqual(r.valid_mean.tolist(),old.valid_mean.tolist())
        self.assertNotIn('tide_lag_+2.0h',r)
        self.assertEqual(r.tide_height_m.tolist(),[3.,2.])
        new.loc[1,'tide_height_m']=np.nan
        with self.assertRaises(ValueError):refresh(p,new)
        r,a=refresh(p,new,True)
        self.assertEqual(r.date.tolist(),['2024-01-01']);self.assertFalse(a.matched.iloc[1])
    def test_groups_month_balance_ties_and_determinism(self):
        d=pd.DataFrame({'date':['2024-01-01','2024-01-02','2024-01-03','2024-02-01','2024-02-02','2024-02-03'],
                        'tide_height_m':[1,2,3,4,5,6]})
        hi,lo=select_groups(d)
        self.assertEqual(hi,['2024-01-03','2024-02-03']);self.assertEqual(lo,['2024-01-01','2024-02-01'])
        self.assertEqual((hi,lo),select_groups(d.sample(frac=1,random_state=1)))
        d.tide_height_m=2
        with self.assertRaises(ValueError):select_groups(d)
    def test_js_only_requested_array_changes(self):
        t="var DATES=['2024-01-01'];\nvar LOW_DATES=['2024-01-01'];\nvar X=42;"
        out,report=replace_array(t,'LOW_DATES',['2024-01-02'])
        self.assertIn("var DATES=['2024-01-01'];",out)
        self.assertIn('var X=42;',out)
        self.assertTrue(report['changed'])
    def test_waterline_equal_levels_cannot_be_split(self):
        tide=np.array([1,1,2,2.]);wet=np.array([0,1,0,1])
        h,kind,err=waterline_threshold(tide,wet)
        self.assertAlmostEqual(err,np.mean((tide>h)!=wet) if kind=='resolved' else .5)
        # A constant datum translation translates the inferred level and leaves error unchanged.
        shifted=waterline_threshold(tide+3,wet)
        self.assertAlmostEqual(shifted[0],h+3);self.assertEqual(shifted[1:],(kind,err))
    def test_bad_date_tables(self):
        p=self.root/'date.csv'
        pd.DataFrame({'date':['2024-01-01']*2,'tide_height_m':[1,2]}).to_csv(p,index=False)
        with self.assertRaises(ValueError):read_date_table(p)

class CliIntegrationTests(unittest.TestCase):
    setUp = TideMigrationTests.setUp
    tearDown = TideMigrationTests.tearDown
    def run_script(self, name, *args):
        return subprocess.run([sys.executable, str(Path(__file__).resolve().parents[1]/name), *map(str,args)],
                              capture_output=True,text=True,check=True)
    def test_required_tiles_and_strict_dates(self):
        sc=self.root/'scenes.csv';td=self.root/'tides.csv';out=self.root/'prepared'
        pd.DataFrame({'scene_id':['a','b'],'overpass_utc':['2024-01-01T03:00Z','2024-01-02T03:00Z'],
                      'cloud_pct':[1,2],'mangrove_valid_frac':[.9,.7],'platform':['S2A','S2A'],
                      'tile_id':['tileA','tileA'],'solar_zenith':[30,32]}).to_csv(sc,index=False)
        pd.DataFrame({'overpass_utc':['2024-01-01T03:00Z','2024-01-02T03:00Z'],
                      'tide_height_m':[2.,3.],'tide_source':'observed_station','tide_datum':'station zero'}).to_csv(td,index=False)
        self.run_script('prepare_dates.py','--scene-list',sc,'--tides',td,'--out-dir',out,
                        '--required-tiles','1','--min-valid','.6','--strict-min-valid','.8')
        d=pd.read_csv(out/'analysis_dates_all.csv')
        self.assertEqual(len(d),2)
        self.assertEqual(d.qualified_strict.tolist(),[True,False])
    def test_cangio_entry_discovers_original_date_names(self):
        site=self.root/'site';(site/'inputs').mkdir(parents=True);(site/'dates').mkdir()
        pd.DataFrame({'scene_id':['a','b','c'],'overpass_utc':['2024-01-01T03:30Z','2024-01-02T03:30Z','2024-01-03T03:30Z']}).to_csv(site/'inputs/cg_scene_list_fullyear.csv',index=False)
        (site/'inputs/h383 (1).csv').write_text('2024,1,1,3,1000\n2024,1,1,4,2000\n2024,1,2,3,2000\n2024,1,2,4,3000\n2024,1,3,3,3000\n2024,1,3,4,4000\n')
        d=pd.DataFrame({'date':['2024-01-01','2024-01-02','2024-01-03'],'tide_height_m':[-1,0,1],
                        'qualified_strict':True,'valid_mean':.9,'month':1})
        d.to_csv(site/'dates/cg_analysis_dates_relaxed06_corr.csv',index=False)
        d.to_csv(site/'dates/cg_dates_all_quality.csv',index=False)
        self.run_script('run_cangio_observed.py','--site-dir',site,'--gee-dir',Path(__file__).resolve().parents[1])
        report=json.loads((site/'observed/migration_inputs.json').read_text())
        self.assertTrue(report['original_analysis_dates'].endswith('cg_analysis_dates_relaxed06_corr.csv'))
        updated=pd.read_csv(site/'observed/dates/cg_analysis_dates_relaxed06_corr.csv')
        self.assertEqual(updated.tide_height_m.tolist(),[1.5,2.5,3.5])
        self.assertTrue(updated.qualified_strict.all())

    def test_sync_report_changed_then_unchanged(self):
        p=self.root/'dates.csv'
        pd.DataFrame({'date':['2024-01-01','2024-01-02','2024-01-03'],
                      'tide_height_m':[1,2,3],'tide_source':'observed_station','tide_datum':'station zero'}).to_csv(p,index=False)
        a=self.root/'first';b=self.root/'second'
        self.run_script('sync_gee_tide_dates.py','--dates',p,'--gee-dir',Path(__file__).resolve().parents[1],'--out-dir',a)
        self.run_script('sync_gee_tide_dates.py','--dates',p,'--gee-dir',a,'--out-dir',b)
        first=json.loads((a/'rerun_plan.json').read_text())['rerun']
        second=json.loads((b/'rerun_plan.json').read_text())['rerun']
        self.assertTrue(first['g08_rerun'])
        self.assertTrue(first['g09_RF_retrain_and_export_all_selected_dates'])
        self.assertFalse(second['g08_rerun'])
        self.assertFalse(second['g09_RF_retrain_and_export_all_selected_dates'])
        self.assertIn("var TIDE_DATE_SOURCE = 'observed_station';",(b/'cg_g08_tide_composite_experiment.js').read_text())

if __name__=='__main__':unittest.main()
