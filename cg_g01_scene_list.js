/**** s01–s02 全年版：越南芹椰影像编目（GEE JavaScript）
 * 与 s01_s02_scene_list_gee_v2.js 完全相同，只把月份扩展到1–12月，导出文件名改为 scene_list_fullyear。
 * 用途：①统计全年每月合格日期数；②为全年面板分析匹配FES潮位（s03_tide_match.py）。
 ****/
// ============ 参数（与项目约定一致，不要改动 bbox 数值）============
var bbox = ee.Geometry.Rectangle([106.43, 10.28, 107.36, 10.94]); // 越南芹椰
var seasonMonths = [1,2,3,4,5,6,7,8,9,10,11,12];  // 全年
var yearStart = 2019, yearEnd = 2024;
// AOI位于4个瓦片交界处，不能先按整幅瓦片云量删瓦片，否则会造成当天空间覆盖不完整。
// 保留全部同日瓦片，后续依据红树林范围内的逐像元有效比例筛日期。
var maxSceneCloud = 100;
var gmwIC = 'projects/sat-io/open-datasets/GMW/extent/GMW_V3';
// ================================================

var aoi = bbox;
Map.centerObject(aoi, 10);

// --- 核验①：打印 AOI 四至和面积，肉眼确认是越南芹椰而不是珠江口 ---
print('=== AOI 核验 ===');
print('AOI 边界坐标 (经纬度)：', aoi.bounds());
print('AOI 面积 (km^2)：', aoi.area().divide(1e6));
print('参考：越南芹椰范围应约在[106.43, 10.28, 107.36, 10.94] 范围内，'
      + '面积约十几到二十几平方公里；若面积达到几千平方公里，说明 AOI 被误改成了大范围（比如珠江口）。');

// --- 月份筛 ---
var monthFilter = ee.Filter.or.apply(null, seasonMonths.map(function (m) {
  return ee.Filter.calendarRange(m, m, 'month');
}));

// --- 红树林掩膜（GMW，仅用于估算每景红树林区有效像元比例）---
// 掩膜：所有年份并集(.max)——不依赖年份属性，覆盖最全。
// 注意：不要改用 .mosaic()，会漏红树林（历史踩过的坑，见 HANDOVER 文档）。
var extent = ee.ImageCollection(gmwIC).max().gt(0).unmask(0).clip(aoi);
var mangrove = extent.unmask(0).gt(0).selfMask();

// --- 核验②：打印红树林掩膜面积，与当地公布的红树林面积对照量级 ---
var mangroveAreaHa = mangrove.multiply(ee.Image.pixelArea()).reduceRegion({
  reducer: ee.Reducer.sum(), geometry: aoi, scale: 30, maxPixels: 1e13, tileScale: 8
}).values().get(0);
print('红树林掩膜面积 (ha)（与当地公布面积对照量级）：', ee.Number(mangroveAreaHa).divide(1e4));

// --- 影像集 ---
var coll = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi)
  .filter(ee.Filter.calendarRange(yearStart, yearEnd, 'year'))
  .filter(monthFilter)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', maxSceneCloud));

print('=== 影像集核验 ===');
print('候选影像总数：', coll.size());

// 核验③：按年份统计影像数，确认 2019-2024 都有覆盖，没有异常集中在某一年
var yearCounts = ee.List.sequence(yearStart, yearEnd).map(function (y) {
  y = ee.Number(y);
  var n = coll.filter(ee.Filter.calendarRange(y, y, 'year')).size();
  return ee.String(y.format()).cat(': ').cat(n.format());
});
print('按年份影像数：', yearCounts);

// 核验④：按瓦片(MGRS tile)统计影像数，确认只涉及越南芹椰附近少数瓦片
var tileList = coll.aggregate_array('MGRS_TILE').distinct();
print('涉及的 MGRS 瓦片：', tileList);
print('参考：记下研究区涉及的瓦片数（pipeline.py dates 会自动按“当天瓦片齐全”筛选）；'
      + '如果出现其它瓦片(如更多远处瓦片)，说明 AOI 比预期大，需要检查 bbox 是否被改动。');

// --- SCL 云掩膜，用于估每景红树林区的有效(非云)像元比例 ---
function sclGood(img) {
  var scl = img.select('SCL');
  return scl.eq(4).or(scl.eq(5)).or(scl.eq(6));
}

// --- 逐景提取属性：ID、过境时刻、云量、红树林区有效像元比例 ---
function sceneInfo(img) {
  img = ee.Image(img);
  var t = ee.Date(img.get('system:time_start'));
  var overpassUtc = t.format('YYYY-MM-dd HH:mm:ss');   // UTC，不做时区转换
  var good = sclGood(img).updateMask(mangrove);
  var validFrac = good.reduceRegion({
    reducer: ee.Reducer.mean(), geometry: aoi, scale: 60,
    maxPixels: 1e13, tileScale: 8
  }).get('SCL');
  return ee.Feature(null, {
    scene_id: img.get('system:index'),
    overpass_utc: overpassUtc,
    cloud_pct: img.get('CLOUDY_PIXEL_PERCENTAGE'),
    mangrove_valid_frac: validFrac,
    platform: img.get('SPACECRAFT_NAME'),
    tile_id: img.get('MGRS_TILE'),
    processing_baseline: img.get('PROCESSING_BASELINE'),
    solar_zenith: img.get('MEAN_SOLAR_ZENITH_ANGLE')
  });
}

var fc = ee.FeatureCollection(coll.map(sceneInfo)).sort('overpass_utc');

print('=== 影像清单预览（前50行）===');
print(fc.limit(50));

// --- 导出 CSV 到 Google Drive ---
// 用新文件名，不覆盖旧 scene_list.csv，方便新旧对比核验
Export.table.toDrive({
  collection: fc,
  description: 'cg_scene_list_fullyear',      // 导出 scene_list_fullyear.csv
  fileFormat: 'CSV',
  selectors: ['scene_id', 'overpass_utc', 'cloud_pct', 'mangrove_valid_frac',
              'platform', 'tile_id', 'processing_baseline', 'solar_zenith']
});

Map.addLayer(mangrove, { palette: ['1a9641'] }, '红树林掩膜');
Map.addLayer(aoi, {}, 'AOI 越南芹椰', false);
