// 内置日期为旧FES结果；请运行sync_gee_tide_dates.py生成实测版本后再执行。
var TIDE_DATE_SOURCE = 'legacy_fes';
if (TIDE_DATE_SOURCE !== 'observed_station') { throw new Error('先运行sync_gee_tide_dates.py更新实测高低潮日期分组'); }
/**** s16：两期变化检测模拟——逐景红树林制图面积随过境潮位的变化（GEE JavaScript）
 * 站点：越南芹椰（由 make_s16_scripts.py 从 code/s16_template.js 生成，请勿手改）
 *
 * 做法：对每个合格日期，用当天的Sentinel-2影像（同日多瓦片取中位数，SCL去云）单独出红树林图，
 *      按v4林缘类型统计“有效像元数”和“判为红树林的像元数”。本地 s16_change_simulation.py 再换算成面积，
 *      估计潮位造成的面积偏差，并模拟“前期 vs 后期”两期变化检测中的虚假变化。
 * 制图方法（与s15一致）：
 *   R1：NDVI≥0.4 且 MNDWI<0（含水体掩膜）；R2：仅NDVI≥0.4；
 *   RF：随机森林分类，用低潮合成影像训练（v4红树林为标签），逐景应用——对应“一个模型用于多期影像”。
 * 统计范围：v4红树林（0–4类）＋林外300 m评估区（255）。
 * 输出（Drive）：cg_s16_daily_maps.csv
 ****/
var ASSET_ROOT = 'projects/mangrove-494808/assets/';
var PREFIX = 'cg_';
var aoi = ee.Geometry.Rectangle([106.43, 10.28, 107.36, 10.94]);
var PROJ = ee.Projection('EPSG:32648').atScale(10);
var ZONE_M = 300;
var DATES = [
 '2019-01-13','2019-01-18','2019-02-07','2019-02-12','2019-02-17','2019-02-27','2019-03-19','2019-04-08',
 '2019-04-23','2019-07-07','2019-07-12','2019-07-22','2019-07-27','2019-08-31','2019-09-25','2019-10-05',
 '2019-10-15','2019-10-25','2019-10-30','2019-11-19','2019-11-24','2019-11-29','2019-12-09','2019-12-14',
 '2019-12-24','2019-12-29','2020-01-03','2020-01-13','2020-01-18','2020-01-23','2020-01-28','2020-02-02',
 '2020-02-12','2020-02-22','2020-02-27','2020-04-02','2020-05-02','2020-06-01','2020-08-25','2020-11-18',
 '2020-11-23','2020-12-03','2021-01-07','2021-01-12','2021-01-17','2021-02-01','2021-02-16','2021-03-03',
 '2021-03-08','2021-03-13','2021-03-23','2021-03-28','2021-04-07','2021-04-12','2021-05-02','2021-06-11',
 '2021-07-01','2021-11-08','2021-12-08','2021-12-23','2022-01-07','2022-01-22','2022-02-16','2022-02-21',
 '2022-02-26','2022-03-08','2022-04-22','2022-06-06','2022-06-11','2022-12-08','2022-12-18','2022-12-23',
 '2022-12-28','2023-02-06','2023-02-11','2023-02-16','2023-02-21','2023-02-26','2023-03-08','2023-03-13',
 '2023-03-18','2023-04-02','2023-05-02','2023-05-22','2023-06-11','2023-10-04','2024-01-12','2024-01-17',
 '2024-01-22','2024-01-27','2024-02-01','2024-02-06','2024-02-11','2024-02-16','2024-02-26','2024-03-02',
 '2024-03-22','2024-04-11','2024-04-16','2024-04-26','2024-05-01','2024-08-09','2024-10-13','2024-11-17',
 '2024-11-22'];
// 越南芹椰：同月配对高潮/低潮日期（每月上/下三分之一），由 make_tide_composite_dates.py 生成
// 高潮组均值 +0.23 m，低潮组均值 -1.46 m（相对平均海平面）
var LOW_DATES = [
 '2019-02-07','2019-04-08','2019-07-22','2019-10-30','2019-12-14','2019-12-29','2020-01-13','2020-01-28',
 '2020-02-12','2020-08-25','2020-11-18','2021-01-17','2021-02-01','2021-03-03','2021-03-13','2021-03-28',
 '2021-04-12','2021-06-11','2021-11-08','2021-12-08','2023-02-06','2023-02-21','2023-03-08','2023-05-22',
 '2023-10-04','2024-01-12','2024-01-27','2024-02-11','2024-02-26','2024-04-11'];

Map.setCenter(106.895, 10.61, 13);
print('逐景日期数', DATES.length, '训练用低潮日期数', LOW_DATES.length);

var CLASS_BANDS = ['edge_class', 'd_boundary', 'd_tidal', 'd_nontidal', 'outside_class', 'gmw_edge_zone'];
var cls = ee.Image(ASSET_ROOT + PREFIX + 'edge_class_v4').rename(CLASS_BANDS);
var edgeAll = cls.select('edge_class').round().toInt();
var isM = edgeAll.neq(255).rename('isM');
var dM = isM.fastDistanceTransform(64, 'pixels', 'squared_euclidean').sqrt().multiply(10).unmask(9999).reproject(PROJ);
var zone = dM.lte(ZONE_M).clip(aoi.buffer(-100)).selfMask();
var grp = edgeAll.rename('grp');

var BANDS = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12'];
var FEATS = BANDS.concat(['NDVI', 'EVI', 'NDRE', 'MNDWI', 'NIRv']);
function maskScl(img) { var s = img.select('SCL'); return img.updateMask(s.eq(4).or(s.eq(5)).or(s.eq(6))); }
function addIdx(x) {
  var ndvi = x.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var evi = x.expression('2.5*(N-R)/(N+6*R-7.5*B+1)', {N: x.select('B8'), R: x.select('B4'), B: x.select('B2')}).rename('EVI');
  return x.addBands([ndvi, evi, x.normalizedDifference(['B6', 'B5']).rename('NDRE'),
    x.normalizedDifference(['B3', 'B11']).rename('MNDWI'), ndvi.multiply(x.select('B8')).rename('NIRv')]);
}
var S2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(aoi);
function dayFilter(dates) {
  return ee.Filter.or.apply(null, dates.map(function (s) { var d = ee.Date(s); return ee.Filter.date(d, d.advance(1, 'day')); }));
}
// RF分类器：低潮合成影像训练
var lowComp = addIdx(S2.filter(dayFilter(LOW_DATES)).map(maskScl).select(BANDS).median().divide(1e4));
var labPts = isM.rename('lab').updateMask(zone).stratifiedSample({numPoints: 1500, classBand: 'lab', region: aoi,
  projection: PROJ, scale: 10, seed: 16, geometries: true, tileScale: 16});
var rfc = ee.Classifier.smileRandomForest({numberOfTrees: 100, seed: 16})
  .train(lowComp.select(FEATS).sampleRegions({collection: labPts, properties: ['lab'], projection: PROJ, scale: 10, tileScale: 16}), 'lab', FEATS);

function daily(ds) {
  ds = ee.String(ds);
  var d0 = ee.Date(ds);
  var x = addIdx(S2.filterDate(d0, d0.advance(1, 'day')).map(maskScl).select(BANDS).median().divide(1e4));
  var valid = x.select('B4').mask().gt(0).rename('n_valid');
  var r1 = x.select('NDVI').gte(0.4).and(x.select('MNDWI').lt(0)).rename('n_R1');
  var r2 = x.select('NDVI').gte(0.4).rename('n_R2');
  var rf = x.select(FEATS).classify(rfc).eq(1).rename('n_RF');
  var stack = valid.addBands([r1.unmask(0), r2.unmask(0), rf.unmask(0)]).updateMask(valid).addBands(grp).updateMask(zone);
  var r = stack.reduceRegion({reducer: ee.Reducer.sum().repeat(4).group({groupField: 4, groupName: 'grp'}),
    geometry: aoi, crs: PROJ, maxPixels: 1e13, tileScale: 16});
  return ee.List(r.get('groups')).map(function (g) {
    g = ee.Dictionary(g); var s = ee.List(g.get('sum'));
    return ee.Feature(null, {date: ds, grp: g.get('grp'), n_valid: s.get(0), n_R1: s.get(1), n_R2: s.get(2), n_RF: s.get(3)});
  });
}
var out = ee.FeatureCollection(ee.List(DATES).map(daily).flatten());
// v4各类面积（像元数），本地换算面积用
var total = isM.rename('n_ref').addBands(grp).updateMask(zone).reduceRegion({
  reducer: ee.Reducer.count().group({groupField: 1, groupName: 'grp'}), geometry: aoi, crs: PROJ, maxPixels: 1e13, tileScale: 16});
var refRows = ee.FeatureCollection(ee.List(total.get('groups')).map(function (g) {
  g = ee.Dictionary(g);
  return ee.Feature(null, {date: 'V4_TOTAL', grp: g.get('grp'), n_valid: g.get('count'), n_R1: -1, n_R2: -1, n_RF: -1});
}));
Export.table.toDrive({collection: out.merge(refRows), description: PREFIX + 's16_daily_maps', fileFormat: 'CSV',
  selectors: ['date', 'grp', 'n_valid', 'n_R1', 'n_R2', 'n_RF']});
