/**** s07 v4：按v4林缘类型固定点提取逐日光谱（GEE JavaScript）
 * 输入：s06 v4 导出的点资产 cg_edge_points_v4；SELECTED_DATES 由本地
 *       python make_gee_date_list.py --dates <合格日期表> 生成后整段粘贴替换。
 * 默认先放原35个10–12月日期，可直接与v3结果对照；全年日期表生成后再替换重跑。
 * 输出列与v3一致，另加 edge_class、d_boundary、d_tidal、d_nontidal、gmw_edge_zone 等。
 ****/
var aoi = ee.Geometry.Rectangle([106.43, 10.28, 107.36, 10.94]);
var ASSET_ROOT = 'projects/mangrove-494808/assets/'; // GEE资产目录（站点配置 asset_root）
var OUT_NAME = 'cg_edge_spectra_v4_fullyear';      // 换全年日期时改名，如 cg_edge_spectra_v4_fullyear
var SELECTED_DATES = [
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
if (SELECTED_DATES.length === 0) { throw new Error('日期尚未写入：先完成编目→潮位→日期步骤，再运行 python pipeline.py gee 重新生成本脚本'); }


var samples = ee.FeatureCollection(ASSET_ROOT + 'cg_edge_points_v4');
Map.setCenter(106.895, 10.61, 13);  // 不访问服务器，避免网络不稳时报错
Map.addLayer(samples, {color: 'cyan'}, 'v4采样点');
print('点数', samples.size(), '各类型', samples.aggregate_histogram('edge_class'));

function maskScl(img) {
  var s = img.select('SCL');
  return img.updateMask(s.eq(4).or(s.eq(5)).or(s.eq(6)));
}
var bands = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12'];
function dayImage(dateStr) {
  var d0 = ee.Date(dateStr), d1 = d0.advance(1, 'day');
  var raw = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(aoi).filterDate(d0, d1);
  var x = raw.map(maskScl).median();
  var nir = x.select('B8').divide(1e4), red = x.select('B4').divide(1e4), blue = x.select('B2').divide(1e4);
  var ndvi = x.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var evi = nir.subtract(red).multiply(2.5)
    .divide(nir.add(red.multiply(6)).subtract(blue.multiply(7.5)).add(1)).rename('EVI');
  var ndre = x.normalizedDifference(['B6', 'B5']).rename('NDRE');
  var mndwi = x.normalizedDifference(['B3', 'B11']).rename('MNDWI');
  return x.select(bands).addBands([ndvi, evi, ndre, mndwi]).set({
    n_source_tiles: raw.aggregate_count_distinct('MGRS_TILE'),
    platform: raw.first().get('SPACECRAFT_NAME'),
    cloud_mean: raw.aggregate_mean('CLOUDY_PIXEL_PERCENTAGE'),
    solar_zenith: raw.aggregate_mean('MEAN_SOLAR_ZENITH_ANGLE')});
}
var out = ee.FeatureCollection(ee.List(SELECTED_DATES).map(function (ds) {
  ds = ee.String(ds);
  var im = dayImage(ds);
  return im.sampleRegions({collection: samples, scale: 10, projection: 'EPSG:32648',
    geometries: false, tileScale: 8})
    .map(function (f) { return f.set({date: ds, n_source_tiles: im.get('n_source_tiles'),
      platform: im.get('platform'), cloud_mean: im.get('cloud_mean'), solar_zenith: im.get('solar_zenith')}); });
})).flatten();

print('日期数', SELECTED_DATES.length, '输出行数', out.size());
print('预览', out.limit(5));
Export.table.toDrive({collection: out, description: OUT_NAME, fileFormat: 'CSV',
  selectors: ['pt_id', 'lon', 'lat', 'date', 'n_source_tiles', 'platform', 'cloud_mean', 'solar_zenith',
    'stratum', 'edge_class', 'shore_id', 'block_id', 'd_boundary', 'd_tidal', 'd_nontidal', 'gmw_edge_zone',
    'optical_water_freq', 'valid_obs', 'water_freq_allyear', 'ndvi_median_allyear', 'canopy_ht',
    'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12', 'NDVI', 'EVI', 'NDRE', 'MNDWI']});
