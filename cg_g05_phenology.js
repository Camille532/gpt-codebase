/**** s10：按v4林缘类型统计全年逐月NDVI物候（GEE JavaScript）
 * 目的：①看全年各月NDVI差异有多大（亚热带红树林是否可以不限定10–12月）；
 *       ②看核心/过渡/向海/向陆林缘的季节曲线是否不同——若不同，全年面板必须加“林缘类型×月份”控制。
 * 输入：s06_edge_classify_v4_gee.js 导出的 cg_edge_class_v4 栅格资产。
 * 输出：两张表（Drive）
 *   cg_phenology_month_by_class      ：多年合并的逐月中位数合成 → 各类型NDVI均值/四分位
 *   cg_phenology_yearmonth_by_class  ：逐年逐月中位数合成 → 各类型NDVI均值（看年际波动、算季节振幅的不确定性）
 ****/
var aoi = ee.Geometry.Rectangle([106.43, 10.28, 107.36, 10.94]);
var yearStart = 2019, yearEnd = 2024;
var ASSET_ROOT = 'projects/mangrove-494808/assets/'; // GEE资产目录（站点配置 asset_root）
var PROJ = ee.Projection('EPSG:32648').atScale(10);
var CLASS_NAMES = ['core', 'transition', 'seaward_edge', 'landward_edge', 'mixed_edge'];

// 按波段顺序重命名（兼容上传的GeoTIFF）；edge_class=255表示非红树林
var CLASS_BANDS = ['edge_class', 'd_boundary', 'd_tidal', 'd_nontidal', 'outside_class', 'gmw_edge_zone'];
var cls = ee.Image(ASSET_ROOT + 'cg_edge_class_v4').rename(CLASS_BANDS);
cls = cls.addBands(cls.select('edge_class').updateMask(cls.select('edge_class').neq(255)), null, true);
cls = cls.select('edge_class');

function maskScl(img) {
  var s = img.select('SCL');
  return img.updateMask(s.eq(4).or(s.eq(5)).or(s.eq(6)));
}
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(aoi)
  .filter(ee.Filter.calendarRange(yearStart, yearEnd, 'year')).map(maskScl)
  .map(function (img) { return img.normalizedDifference(['B8', 'B4']).rename('NDVI')
    .copyProperties(img, ['system:time_start']); });

// 按林缘类型分组统计；中位数合成可压低单日高潮、薄云等异常
function statsByClass(ndvi, props) {
  var grouped = ndvi.addBands(cls).reduceRegion({
    reducer: ee.Reducer.mean().combine(ee.Reducer.percentile([25, 75]), '', true)
      .combine(ee.Reducer.count(), '', true).group({groupField: 1, groupName: 'edge_class'}),
    geometry: aoi, crs: PROJ, maxPixels: 1e13, tileScale: 16});
  return ee.List(grouped.get('groups')).map(function (g) {
    g = ee.Dictionary(g);
    return ee.Feature(null, g.combine(props)
      .set('class_name', ee.List(CLASS_NAMES).get(g.getNumber('edge_class').int())));
  });
}

var months = ee.List.sequence(1, 12);
var monthRows = months.map(function (m) {
  var c = s2.filter(ee.Filter.calendarRange(m, m, 'month'));
  return statsByClass(c.median(), ee.Dictionary({month: m, n_images: c.size()}));
}).flatten();
var monthFc = ee.FeatureCollection(monthRows);

var ymRows = ee.List.sequence(yearStart, yearEnd).map(function (y) {
  return months.map(function (m) {
    var c = s2.filter(ee.Filter.calendarRange(y, y, 'year')).filter(ee.Filter.calendarRange(m, m, 'month'));
    // 某年某月无影像时返回空列表
    return ee.Algorithms.If(c.size().gt(0),
      statsByClass(c.median(), ee.Dictionary({year: y, month: m, n_images: c.size()})), ee.List([]));
  }).flatten();
}).flatten();
var ymFc = ee.FeatureCollection(ymRows);

print('逐月×林缘类型', monthFc);
print(ui.Chart.feature.groups(monthFc.filter(ee.Filter.lte('edge_class', 3)), 'month', 'mean', 'class_name')
  .setChartType('LineChart').setOptions({title: '越南芹椰红树林逐月NDVI（按v4林缘类型，2019–2024合并）',
    hAxis: {title: '月份', ticks: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]}, vAxis: {title: 'NDVI均值'},
    pointSize: 4}));

Export.table.toDrive({collection: monthFc, description: 'cg_phenology_month_by_class',
  fileFormat: 'CSV', selectors: ['month', 'edge_class', 'class_name', 'mean', 'p25', 'p75', 'count', 'n_images']});
Export.table.toDrive({collection: ymFc, description: 'cg_phenology_yearmonth_by_class',
  fileFormat: 'CSV', selectors: ['year', 'month', 'edge_class', 'class_name', 'mean', 'p25', 'p75', 'count', 'n_images']});
