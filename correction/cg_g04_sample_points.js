/**** s06c v4（第3步/共3步）：检查分类结果并分层采样（GEE JavaScript）
 * 输入：s06a 的 cg_base_v4、s06b 的 cg_edge_class_v4（都已是资产，计算很轻）。
 * 分层：edge_class(0核心 1过渡 2向海 3向陆；4混合不采样) × 岸别(北/南)，每层PTS_PER_STRATUM点。
 * 操作：Run → 看Console面积和各类型统计 → Tasks运行两个导出：
 *   cg_edge_sample_points_v4（CSV，放 tidedata\cg_edge_v4_inputs\）
 *   cg_edge_points_v4_asset（资产，s07 v4读取）
 ****/
var aoi = ee.Geometry.Rectangle([106.43, 10.28, 107.36, 10.94]);
var gmwIC = 'projects/sat-io/open-datasets/GMW/extent/GMW_V3';
var PROJ = ee.Projection('EPSG:32648').atScale(10);
var ASSET_ROOT = 'projects/mangrove-494808/assets/'; // GEE资产目录（站点配置 asset_root）
var PTS_PER_STRATUM = 150, SEED = 2026, MIN_VALID_OBS = 20;
var shoreSplit = 106.83, blockSizeDeg = 0.01;

// 按波段顺序重命名：直接导出的资产与“Drive导出GeoTIFF再上传”的资产（波段名b1..b7）都能读取；-9999为无数据
var BASE_BANDS = ['ndvi_med', 'mndwi_med', 'b11_med', 'water_freq_all', 'valid_all', 'water_freq_ond', 'valid_ond'];
var base = ee.Image(ASSET_ROOT + 'cg_base_v4').rename(BASE_BANDS);
base = base.updateMask(base.neq(-9999));
// 按波段顺序重命名（兼容上传的GeoTIFF）；edge_class=255表示非红树林
var CLASS_BANDS = ['edge_class', 'd_boundary', 'd_tidal', 'd_nontidal', 'outside_class', 'gmw_edge_zone'];
var cls = ee.Image(ASSET_ROOT + 'cg_edge_class_v4').rename(CLASS_BANDS);
cls = cls.addBands(cls.select('edge_class').updateMask(cls.select('edge_class').neq(255)), null, true);
var edgeClass = cls.select('edge_class').round().toInt();
var mangrove = edgeClass.mask().gt(0).unmask(0);
var gmw = ee.ImageCollection(gmwIC).max().gt(0).unmask(0).clip(aoi).reproject(PROJ);

// ---------- 面积检查 ----------
function areaHa(mask) {
  return ee.Number(mask.selfMask().multiply(ee.Image.pixelArea()).reduceRegion({
    reducer: ee.Reducer.sum(), geometry: aoi, crs: PROJ, maxPixels: 1e13, tileScale: 8}).values().get(0)).divide(1e4);
}
print('GMW并集面积 ha', areaHa(gmw));
print('红树林v4面积 ha', areaHa(mangrove));
print('GMW内被剔除 ha（泥滩/水面等）', areaHa(gmw.and(mangrove.not())));
print('GMW外新增 ha', areaHa(mangrove.and(gmw.not())));
print('潮汐连通区面积 ha', areaHa(cls.select('outside_class').eq(1)));
print('各林缘类型像元数（0核心 1过渡 2向海 3向陆 4混合）', edgeClass.reduceRegion({
  reducer: ee.Reducer.frequencyHistogram(), geometry: aoi, crs: PROJ, maxPixels: 1e13, tileScale: 8}));

// ---------- 属性与采样 ----------
var validOnd = base.select('valid_ond').rename('valid_obs');
var canopy = ee.ImageCollection('projects/sat-io/open-datasets/global_mangrove_distribution/hba95')
  .mosaic().rename('canopy_ht').unmask(-9999);
var ll = ee.Image.pixelLonLat(), lon = ll.select('longitude'), lat = ll.select('latitude');
var shore = ll.select('longitude').gt(shoreSplit).rename('shore_id');
var blockId = lat.subtract(10.28).divide(blockSizeDeg).floor().multiply(1000)
  .add(lon.subtract(106.43).divide(blockSizeDeg).floor()).rename('block_id').toInt();
// AOI矩形外被当作非红树林，会在矩形边上形成伪林缘；距AOI边界100 m内不采样
var strata = edgeClass.multiply(10).add(shore).rename('stratum')
  .updateMask(edgeClass.lte(3)).updateMask(validOnd.gte(MIN_VALID_OBS))
  .clip(aoi.buffer(-100)).toInt();

var stack = strata.addBands([edgeClass.rename('edge_class'), shore, blockId,
  cls.select(['d_boundary', 'd_tidal', 'd_nontidal']), cls.select('gmw_edge_zone').round().toInt(),
  base.select('water_freq_ond').rename('optical_water_freq'), validOnd,
  base.select('water_freq_all').rename('water_freq_allyear'), base.select('ndvi_med').rename('ndvi_median_allyear'),
  canopy]);
var samples = stack.stratifiedSample({numPoints: PTS_PER_STRATUM, classBand: 'stratum', region: aoi,
  projection: PROJ, scale: 10, geometries: true, seed: SEED, tileScale: 16});
samples = samples.map(function (f) {
  var c = ee.List(f.geometry().coordinates());
  return f.set({pt_id: f.get('system:index'), lon: c.get(0), lat: c.get(1)});
});
print('采样点数', samples.size());
print('分层×岸别（十位=类型，个位=岸别 0南1北）', samples.aggregate_histogram('stratum'));
print('v4类型_v3距离带 交叉', samples.map(function (f) {
  return f.set('cross', ee.String(f.getNumber('edge_class').format('%d')).cat('_')
    .cat(f.getNumber('gmw_edge_zone').format('%d')));
}).aggregate_histogram('cross'));

Map.setCenter(106.895, 10.61, 13);
Map.setOptions('SATELLITE');
Map.addLayer(edgeClass, {min: 0, max: 4, palette: ['1a9641', 'a6d96a', 'd7191c', 'ff8c00', 'ff00ff']},
  '林缘类型 0核心 1过渡 2向海 3向陆 4混合');
Map.addLayer(samples, {color: 'white'}, 'v4采样点');

var cols = ['pt_id', 'lon', 'lat', 'stratum', 'edge_class', 'shore_id', 'block_id',
  'd_boundary', 'd_tidal', 'd_nontidal', 'gmw_edge_zone', 'optical_water_freq', 'valid_obs',
  'water_freq_allyear', 'ndvi_median_allyear', 'canopy_ht'];
Export.table.toDrive({collection: samples, description: 'cg_edge_sample_points_v4',
  fileFormat: 'CSV', selectors: cols});
Export.table.toAsset({collection: samples, description: 'cg_edge_points_v4_asset',
  assetId: ASSET_ROOT + 'cg_edge_points_v4'});
