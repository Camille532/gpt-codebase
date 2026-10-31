/**** s06a v4（第1步/共3步）：导出全年光谱中位数与水体频率底图（GEE JavaScript）
 * 为什么拆开：2019–2024全年约1700景影像的中位数合成很耗内存，与连通、距离计算放在同一次交互请求里
 *            会报“User memory limit exceeded”。先把它在后台导出为资产，后两步直接读取。
 * 操作：Run → Tasks里运行 cg_base_v4_asset（约10–30分钟）。Console无需等待结果。
 * 输出资产 cg_base_v4 的波段：
 *   ndvi_med, mndwi_med, b11_med       全年逐像元中位数
 *   water_freq_all, valid_all          全年光学水体频率（MNDWI>0）与有效观测次数
 *   water_freq_ond, valid_ond          10–12月水体频率与有效观测次数（与v3属性口径一致）
 ****/
var aoi = ee.Geometry.Rectangle([106.43, 10.28, 107.36, 10.94]);
var yearStart = 2019, yearEnd = 2024;
var ASSET_ROOT = 'projects/mangrove-494808/assets/'; // GEE资产目录（站点配置 asset_root）

function prep(img) {
  var s = img.select('SCL');
  var good = s.eq(4).or(s.eq(5)).or(s.eq(6));
  var x = img.updateMask(good);
  var mndwi = x.normalizedDifference(['B3', 'B11']).rename('MNDWI');
  var valid = x.select('B3').mask().rename('VALID');
  return x.normalizedDifference(['B8', 'B4']).rename('NDVI')
    .addBands([mndwi, x.select('B11'), mndwi.gt(0).unmask(0).updateMask(valid).rename('WATER'), valid])
    .copyProperties(img, ['system:time_start']);
}
var all = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(aoi)
  .filter(ee.Filter.calendarRange(yearStart, yearEnd, 'year'))
  .select(['B3', 'B4', 'B8', 'B11', 'SCL']).map(prep);
var ond = all.filter(ee.Filter.calendarRange(10, 12, 'month'));

var med = all.select(['NDVI', 'MNDWI', 'B11']).median().rename(['ndvi_med', 'mndwi_med', 'b11_med']);
var validAll = all.select('VALID').sum();
var validOnd = ond.select('VALID').sum();
var base = med.addBands([
  all.select('WATER').sum().divide(validAll).rename('water_freq_all'), validAll.rename('valid_all'),
  ond.select('WATER').sum().divide(validOnd).rename('water_freq_ond'), validOnd.rename('valid_ond')
]).toFloat().unmask(-9999);  // 无数据写-9999，读取时再掩膜

print('影像数（全年）', all.size());
Map.setCenter(106.895, 10.61, 13);
Map.addLayer(base.select('ndvi_med'), {min: 0, max: 0.9, palette: ['8c510a', 'f6e8c3', '1a9850']}, 'NDVI全年中位数（预览）');

// 路线①：直接导出资产
Export.image.toAsset({image: base, description: 'cg_base_v4_asset', assetId: ASSET_ROOT + 'cg_base_v4',
  region: aoi, crs: 'EPSG:32648', scale: 10, maxPixels: 1e10});
// 路线②（资产导出失败时用）：导出GeoTIFF到Drive → 下载 → Assets面板 NEW → Image upload 上传，资产名填 cg_base_v4
Export.image.toDrive({image: base, description: 'cg_base_v4_drive', fileNamePrefix: 'cg_base_v4',
  region: aoi, crs: 'EPSG:32648', scale: 10, maxPixels: 1e10});
