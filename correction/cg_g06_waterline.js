/**** s11：红树林向海前缘外侧潮滩的逐日干湿序列（GEE JavaScript）
 * 目的：回答“卫星过境时海水到底有没有到达红树林前缘”，并用水边线法估计前缘潮滩高程。
 * 思路：取v4分类中与湾连通的潮汐区像元，按“到最近红树林的距离”分档；
 *       逐日记录这些像元是否为水（MNDWI>0），本地再与实测或模型过境潮位配对（s11_waterline_elevation.py）。
 *       本地高程估计沿用输入潮位的基准；UHSLC站点零点不可直接当作平均海平面。
 * 使用64个10–12月独立日期（不只35个合格日期）：逐像元SCL掩膜，部分云天也能提供有效像元。
 * 输入：s06 v4导出的 cg_edge_class_v4 栅格资产。
 * 输出：①front_points_daily：逐点逐日MNDWI；②front_bins_daily：逐日各距离档的湿像元比例。
 ****/
var aoi = ee.Geometry.Rectangle([106.43, 10.28, 107.36, 10.94]);
var ASSET_ROOT = 'projects/mangrove-494808/assets/'; // GEE资产目录（站点配置 asset_root）
var PROJ = ee.Projection('EPSG:32648').atScale(10);
var PTS_PER_BIN = 200, SEED = 11;
// 距最近红树林像元的距离档（m，像元中心距离；紧邻红树林的像元为10 m）
var BIN_EDGES = [0, 10, 20, 30, 60, 100, 200, 500];
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
if (DATES.length === 0) { throw new Error('日期尚未写入：先完成编目→潮位→日期步骤，再运行 python pipeline.py gee 重新生成本脚本'); }


// 按波段顺序重命名（兼容上传的GeoTIFF）；edge_class=255表示非红树林
var CLASS_BANDS = ['edge_class', 'd_boundary', 'd_tidal', 'd_nontidal', 'outside_class', 'gmw_edge_zone'];
var cls = ee.Image(ASSET_ROOT + 'cg_edge_class_v4').rename(CLASS_BANDS);
cls = cls.addBands(cls.select('edge_class').updateMask(cls.select('edge_class').neq(255)), null, true);
var mangrove = cls.select('edge_class').mask().gt(0).unmask(0).reproject(PROJ);
var tidal = cls.select('outside_class').eq(1).unmask(0).and(mangrove.not()).reproject(PROJ);
var dM = mangrove.fastDistanceTransform(256, 'pixels', 'squared_euclidean').sqrt().multiply(10)
  .unmask(9999).reproject(PROJ).rename('d_mangrove');

var bin = ee.Image(0);
for (var i = 1; i < BIN_EDGES.length; i++) {
  bin = bin.where(dM.gt(BIN_EDGES[i - 1]).and(dM.lte(BIN_EDGES[i])), i);
}
bin = bin.updateMask(tidal).updateMask(bin.gt(0)).rename('front_bin').toInt()
  .clip(aoi.buffer(-100));
var ll = ee.Image.pixelLonLat();
var shore = ll.select('longitude').gt(106.83).rename('shore_id');

var pts = bin.addBands([dM, shore]).stratifiedSample({numPoints: PTS_PER_BIN, classBand: 'front_bin',
  region: aoi, projection: PROJ, scale: 10, geometries: true, seed: SEED, tileScale: 16});
pts = pts.map(function (f) {
  var c = ee.List(f.geometry().coordinates());
  return f.set({pt_id: f.get('system:index'), lon: c.get(0), lat: c.get(1)});
});
print('前缘采样点', pts.size(), pts.aggregate_histogram('front_bin'));

function maskScl(img) {
  var s = img.select('SCL');
  return img.updateMask(s.eq(4).or(s.eq(5)).or(s.eq(6)));
}
function dayMndwi(ds) {
  var d0 = ee.Date(ds);
  var x = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(aoi)
    .filterDate(d0, d0.advance(1, 'day')).map(maskScl).median();
  return x.normalizedDifference(['B3', 'B11']).rename('MNDWI');
}

var pointRows = ee.FeatureCollection(ee.List(DATES).map(function (ds) {
  ds = ee.String(ds);
  return dayMndwi(ds).sampleRegions({collection: pts, projection: PROJ, scale: 10,
    geometries: false, tileScale: 8}).map(function (f) { return f.set('date', ds); });
})).flatten();

// 面积口径：逐日各距离档中有效像元数与湿像元数（不受采样影响，作为主图）
var binRows = ee.FeatureCollection(ee.List(DATES).map(function (ds) {
  ds = ee.String(ds);
  var m = dayMndwi(ds);
  var wet = m.gt(0).rename('wet');
  var stats = wet.addBands(ee.Image(1).updateMask(m.mask()).rename('valid')).addBands(bin)
    .reduceRegion({reducer: ee.Reducer.sum().repeat(2).group({groupField: 2, groupName: 'front_bin'}),
      geometry: aoi, crs: PROJ, maxPixels: 1e13, tileScale: 16});
  return ee.List(stats.get('groups')).map(function (g) {
    g = ee.Dictionary(g); var s = ee.List(g.get('sum'));
    return ee.Feature(null, {date: ds, front_bin: g.get('front_bin'), n_wet: s.get(0), n_valid: s.get(1)});
  });
}).flatten());

Map.setCenter(106.895, 10.61, 13); Map.setOptions('SATELLITE');  // setCenter不访问服务器
Map.addLayer(bin, {min: 1, max: 7, palette: ['b10026', 'e31a1c', 'fc4e2a', 'fd8d3c', 'feb24c', 'fed976', 'ffffb2']},
  '潮滩距红树林距离档');
Map.addLayer(pts, {color: 'cyan'}, '前缘采样点', false);
print('逐日距离档预览', binRows.limit(10));

Export.table.toDrive({collection: pointRows, description: 'cg_front_points_daily', fileFormat: 'CSV',
  selectors: ['pt_id', 'lon', 'lat', 'front_bin', 'd_mangrove', 'shore_id', 'date', 'MNDWI']});
Export.table.toDrive({collection: binRows, description: 'cg_front_bins_daily', fileFormat: 'CSV',
  selectors: ['date', 'front_bin', 'n_wet', 'n_valid']});
