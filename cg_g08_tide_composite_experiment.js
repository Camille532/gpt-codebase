/**** s15：潮位对红树林制图与生物量/碳储量估算的影响——高潮合成 vs 低潮合成（GEE JavaScript）
 * 站点：越南芹椰（由 make_s15_scripts.py 从 code/s15_template.js 生成，请勿手改）
 *
 * 设计
 *  - 高潮合成、低潮合成：同月配对挑选日期（每月上/下三分之一潮位），两组月份构成相同，差异主要来自潮位。
 *  - 参考生物量：GEE社区目录中Simard等全球红树林产品（优先AGB；找不到则用冠层高度hba95）。
 *    参考值只用来让模型形式符合常见做法；我们关心的是“同一流程换成高潮影像后结果变了多少”，不依赖参考值的绝对精度。
 *  - 生物量模型（常见做法）：
 *      RF回归（Sentinel-2各波段+指数）；线性模型 ref~NDVI、ref~EVI、ref~NIRv。
 *    情景：①低潮训练→分别用于低潮/高潮影像（模型固定、影像潮位不同，对应模型移植与变化检测）；
 *          ②高潮训练→用于高潮影像（整个流程都用高潮影像）。
 *  - 制图：规则①NDVI≥0.4且MNDWI<0（常见的水体掩膜流程）；规则②仅NDVI≥0.4；RF分类（v4红树林为标签）。
 *  - 汇总：v4红树林范围内各林缘类型的总生物量（固定范围，只看密度偏差）；各制图结果的面积及其×生物量（面积+密度合计偏差）。
 * 输出（Drive）：cg_s15_summary.csv（各情景汇总）、cg_s15_pixels.csv（像元样本，供作图）
 * 导出在后台运行，Console只打印少量检查信息，避免交互模式内存超限。
 ****/
var ASSET_ROOT = 'projects/mangrove-494808/assets/';
var PREFIX = 'cg_';
var aoi = ee.Geometry.Rectangle([106.43, 10.28, 107.36, 10.94]);
var PROJ = ee.Projection('EPSG:32648').atScale(10);
var ZONE_M = 300;            // 制图评估区：红树林及其外300 m
var CARBON_FRACTION = 0.47;  // 生物量→碳（IPCC默认值；只做比例换算，不影响相对偏差）
// 越南芹椰：同月配对高潮/低潮日期（每月上/下三分之一），由 make_tide_composite_dates.py 生成
// 高潮组均值 +0.23 m，低潮组均值 -1.46 m（相对平均海平面）
var HIGH_DATES = [
 '2019-01-18','2019-02-27','2019-07-12','2019-08-31','2019-10-25','2019-11-24','2019-12-09','2019-12-24',
 '2020-01-03','2020-01-18','2020-02-02','2020-04-02','2020-05-02','2020-06-01','2020-11-23','2021-01-07',
 '2021-03-08','2021-03-23','2022-02-26','2022-04-22','2022-12-18','2023-02-16','2023-02-26','2023-03-18',
 '2024-01-22','2024-02-06','2024-02-16','2024-03-02','2024-04-16','2024-10-13'];
var LOW_DATES = [
 '2019-02-07','2019-04-08','2019-07-22','2019-10-30','2019-12-14','2019-12-29','2020-01-13','2020-01-28',
 '2020-02-12','2020-08-25','2020-11-18','2021-01-17','2021-02-01','2021-03-03','2021-03-13','2021-03-28',
 '2021-04-12','2021-06-11','2021-11-08','2021-12-08','2023-02-06','2023-02-21','2023-03-08','2023-05-22',
 '2023-10-04','2024-01-12','2024-01-27','2024-02-11','2024-02-26','2024-04-11'];

Map.setCenter(106.895, 10.61, 13);
Map.setOptions('SATELLITE');
print('高潮组日期数', HIGH_DATES.length, '低潮组日期数', LOW_DATES.length);

// ---------- 分类底图 ----------
var CLASS_BANDS = ['edge_class', 'd_boundary', 'd_tidal', 'd_nontidal', 'outside_class', 'gmw_edge_zone'];
var cls = ee.Image(ASSET_ROOT + PREFIX + 'edge_class_v4').rename(CLASS_BANDS);
var edgeAll = cls.select('edge_class').round().toInt();          // 255=非红树林
var isM = edgeAll.neq(255).rename('isM');
var dM = isM.fastDistanceTransform(64, 'pixels', 'squared_euclidean').sqrt().multiply(10).unmask(9999).reproject(PROJ);
var zone = dM.lte(ZONE_M).clip(aoi.buffer(-100)).selfMask();
var grp = edgeAll.rename('grp');                                  // 分组：0–4林缘类型，255=林外

// ---------- 合成影像 ----------
var BANDS = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12'];
var FEATS = BANDS.concat(['NDVI', 'EVI', 'NDRE', 'MNDWI', 'NIRv']);
function maskScl(img) { var s = img.select('SCL'); return img.updateMask(s.eq(4).or(s.eq(5)).or(s.eq(6))); }
function dayFilter(dates) {
  return ee.Filter.or.apply(null, dates.map(function (s) { var d = ee.Date(s); return ee.Filter.date(d, d.advance(1, 'day')); }));
}
function composite(dates) {
  var c = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(aoi).filter(dayFilter(dates)).map(maskScl);
  var x = c.select(BANDS).median().divide(1e4);
  var ndvi = x.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var evi = x.expression('2.5*(N-R)/(N+6*R-7.5*B+1)', {N: x.select('B8'), R: x.select('B4'), B: x.select('B2')}).rename('EVI');
  var ndre = x.normalizedDifference(['B6', 'B5']).rename('NDRE');
  var mndwi = x.normalizedDifference(['B3', 'B11']).rename('MNDWI');
  var nirv = ndvi.multiply(x.select('B8')).rename('NIRv');
  return x.addBands([ndvi, evi, ndre, mndwi, nirv]).set('n_images', c.size());
}
var hi = composite(HIGH_DATES), lo = composite(LOW_DATES);
print('合成影像景数（高/低）', hi.get('n_images'), lo.get('n_images'));

// ---------- 参考生物量 ----------
var MGD = 'projects/sat-io/open-datasets/global_mangrove_distribution';
var refId = null, refType = null, refName = 'canopy_height_hba95_m';
try {
  var la = ee.data.listAssets(MGD);
  print('全球红树林产品目录', (la.assets || []).map(function (a) { return a.id + ' (' + a.type + ')'; }));
  (la.assets || []).forEach(function (a) { if (!refId && /agb/i.test(a.id)) { refId = a.id; refType = a.type; } });
} catch (e) { print('列目录失败，改用冠层高度', e); }
var ref;
if (refId) {
  ref = (refType === 'IMAGE_COLLECTION' ? ee.ImageCollection(refId).mosaic() : ee.Image(refId)).select(0);
  refName = 'agb_Mg_per_ha';
} else {
  ref = ee.ImageCollection(MGD + '/hba95').mosaic().select(0);
}
ref = ref.rename('ref').toFloat();
print('参考变量', refName, refId || (MGD + '/hba95'));

// ---------- 生物量模型 ----------
var trainPts = edgeAll.updateMask(isM).rename('edge_class').addBands(ref).updateMask(ref.gt(0))
  .stratifiedSample({numPoints: 300, classBand: 'edge_class', region: aoi, projection: PROJ, scale: 10,
    seed: 15, geometries: true, tileScale: 16});
function samp(img) {
  return img.select(FEATS).addBands(ref).sampleRegions({collection: trainPts, properties: ['edge_class'],
    projection: PROJ, scale: 10, tileScale: 16});
}
var trLo = samp(lo), trHi = samp(hi);
function rf(tr) {
  return ee.Classifier.smileRandomForest({numberOfTrees: 100, seed: 15}).setOutputMode('REGRESSION').train(tr, 'ref', FEATS);
}
var rfLo = rf(trLo), rfHi = rf(trHi);
function linfit(tr, x) { return ee.Dictionary(tr.reduceColumns(ee.Reducer.linearFit(), [x, 'ref'])); }
function linPred(img, fit, x) { return img.select(x).multiply(fit.getNumber('scale')).add(fit.getNumber('offset')).max(0); }
var fN = linfit(trLo, 'NDVI'), fE = linfit(trLo, 'EVI'), fV = linfit(trLo, 'NIRv');

var preds = {
  'REF': ref.max(0),
  'RF_trainLow_imgLow': lo.select(FEATS).classify(rfLo).max(0),
  'RF_trainLow_imgHigh': hi.select(FEATS).classify(rfLo).max(0),
  'RF_trainHigh_imgHigh': hi.select(FEATS).classify(rfHi).max(0),
  'LIN_NDVI_imgLow': linPred(lo, fN, 'NDVI'), 'LIN_NDVI_imgHigh': linPred(hi, fN, 'NDVI'),
  'LIN_EVI_imgLow': linPred(lo, fE, 'EVI'), 'LIN_EVI_imgHigh': linPred(hi, fE, 'EVI'),
  'LIN_NIRv_imgLow': linPred(lo, fV, 'NIRv'), 'LIN_NIRv_imgHigh': linPred(hi, fV, 'NIRv')
};

// ---------- 制图 ----------
var labPts = isM.rename('lab').updateMask(zone).stratifiedSample({numPoints: 1500, classBand: 'lab', region: aoi,
  projection: PROJ, scale: 10, seed: 16, geometries: true, tileScale: 16});
function sampC(img) { return img.select(FEATS).sampleRegions({collection: labPts, properties: ['lab'], projection: PROJ, scale: 10, tileScale: 16}); }
var rfcLo = ee.Classifier.smileRandomForest({numberOfTrees: 100, seed: 16}).train(sampC(lo), 'lab', FEATS);
var rfcHi = ee.Classifier.smileRandomForest({numberOfTrees: 100, seed: 16}).train(sampC(hi), 'lab', FEATS);
function rule1(img) { return img.select('NDVI').gte(0.4).and(img.select('MNDWI').lt(0)); }
function rule2(img) { return img.select('NDVI').gte(0.4); }
var maps = {
  'R1_NDVI_MNDWImask_imgLow': rule1(lo), 'R1_NDVI_MNDWImask_imgHigh': rule1(hi),
  'R2_NDVIonly_imgLow': rule2(lo), 'R2_NDVIonly_imgHigh': rule2(hi),
  'RFC_trainLow_imgLow': lo.select(FEATS).classify(rfcLo), 'RFC_trainLow_imgHigh': hi.select(FEATS).classify(rfcLo),
  'RFC_trainHigh_imgHigh': hi.select(FEATS).classify(rfcHi)
};

// ---------- 汇总 ----------
var haImg = ee.Image.pixelArea().divide(1e4);
function grouped(img, scenario, metric) {   // img：每像元数值（已乘面积）；按grp分组求和
  var r = img.rename('v').addBands(grp).updateMask(zone).reduceRegion({
    reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'grp'}), geometry: aoi, crs: PROJ, maxPixels: 1e13, tileScale: 16});
  return ee.List(r.get('groups')).map(function (g) {
    g = ee.Dictionary(g);
    return ee.Feature(null, {site: '越南芹椰', scenario: scenario, metric: metric, grp: g.get('grp'), value: g.get('sum'), ref_name: refName});
  });
}
var rows = ee.List([]);
Object.keys(preds).forEach(function (k) {   // 固定v4红树林范围内的总生物量（Mg或m·ha）
  rows = rows.cat(grouped(preds[k].multiply(haImg).updateMask(isM), k, 'biomass_fixed_mask'));
});
Object.keys(maps).forEach(function (k) {    // 制图面积（ha）及该制图范围内的RF生物量（低潮训练模型、同潮位影像）
  var m = maps[k].eq(1);
  var p = (k.indexOf('imgHigh') >= 0) ? preds['RF_trainLow_imgHigh'] : preds['RF_trainLow_imgLow'];
  rows = rows.cat(grouped(m.multiply(haImg), k, 'mapped_area_ha'));
  rows = rows.cat(grouped(p.multiply(haImg).updateMask(m), k, 'biomass_in_mapped_area'));
});
rows = rows.cat(grouped(isM.multiply(haImg), 'V4_REFERENCE', 'mapped_area_ha'));
var summary = ee.FeatureCollection(rows.flatten());

// 模型训练精度（低潮训练集上的拟合，仅供参考）
var fitInfo = ee.FeatureCollection([ee.Feature(null, {site: '越南芹椰', scenario: 'LIN_fit', metric: 'NDVI_scale_offset',
  grp: -1, value: ee.String(fN.getNumber('scale').format('%.4f')).cat(';').cat(fN.getNumber('offset').format('%.4f')), ref_name: refName}),
  ee.Feature(null, {site: '越南芹椰', scenario: 'N_train', metric: 'n_train_low_high', grp: -1,
    value: ee.String(trLo.size().format()).cat(';').cat(trHi.size().format()), ref_name: refName})]);
summary = summary.merge(fitInfo);

// 像元样本（作图用）
var pixPts = edgeAll.updateMask(isM).rename('edge_class').stratifiedSample({numPoints: 400, classBand: 'edge_class',
  region: aoi, projection: PROJ, scale: 10, seed: 17, geometries: true, tileScale: 16});
var pixImg = ref.rename('ref')
  .addBands(preds['RF_trainLow_imgLow'].rename('rf_low')).addBands(preds['RF_trainLow_imgHigh'].rename('rf_high'))
  .addBands(lo.select(['NDVI', 'MNDWI', 'B8', 'B11']).rename(['ndvi_low', 'mndwi_low', 'nir_low', 'swir_low']))
  .addBands(hi.select(['NDVI', 'MNDWI', 'B8', 'B11']).rename(['ndvi_high', 'mndwi_high', 'nir_high', 'swir_high']))
  .addBands(cls.select(['d_boundary', 'd_tidal']));
var pix = pixImg.sampleRegions({collection: pixPts, properties: ['edge_class'], projection: PROJ, scale: 10, tileScale: 16});

Map.addLayer(lo.select(['B4', 'B3', 'B2']), {min: 0, max: 0.15}, '低潮合成 真彩色', false);
Map.addLayer(hi.select(['B4', 'B3', 'B2']), {min: 0, max: 0.15}, '高潮合成 真彩色', false);
Map.addLayer(maps['R1_NDVI_MNDWImask_imgLow'].selfMask(), {palette: ['00ff00']}, '规则①制图 低潮', false);
Map.addLayer(maps['R1_NDVI_MNDWImask_imgHigh'].selfMask(), {palette: ['ff0000']}, '规则①制图 高潮', false);

Export.table.toDrive({collection: summary, description: PREFIX + 's15_summary', fileFormat: 'CSV',
  selectors: ['site', 'scenario', 'metric', 'grp', 'value', 'ref_name']});
Export.table.toDrive({collection: pix, description: PREFIX + 's15_pixels', fileFormat: 'CSV'});
