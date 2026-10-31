/**** s06b v4（第2步/共3步）：红树林边界、潮汐连通区与林缘类型分类（GEE JavaScript）
 * 输入：s06a 导出的 cg_base_v4 资产。
 * 输出资产 cg_edge_class_v4（s06c、s10、s11都读取它），波段：
 *   edge_class   0核心(d_boundary>60 m) 1过渡(30–60 m) 2向海林缘 3向陆林缘 4混合林缘（仅红树林像元有值）
 *   d_boundary   到最近非红树林像元的欧氏距离（m）
 *   d_tidal      到最近潮汐连通区像元的距离（m）
 *   d_nontidal   到最近非潮汐（鱼塘/基围/陆地）像元的距离（m）
 *   outside_class 林外类型 1潮汐区 2非潮汐水体 3非潮汐陆地（红树林内为0）
 *   gmw_edge_zone v3的GMW方形核距离带（0/1/2，GMW外为-1），用于新旧对照
 *
 * 方法要点：
 *  1. 红树林边界：GMW并集外扩30 m内，全年NDVI中位数>=NDVI_MIN且水体频率<=0.5的10 m像元；剔除碎斑，填林内小空洞。
 *  2. 潮汐连通区：林外“水体频率>=5%或湿泥滩光谱”的候选像元中，与开阔湾面连通的部分（cumulativeCost）。
 *     基围、堤岸阻断连通 → 鱼塘归为非潮汐。只用全年光谱，不用分析日期潮位，避免循环论证。
 *  3. 欧氏距离：fastDistanceTransform。
 *
 * 操作：Run → 先看地图图层（放大到局部看，整幅看可能因计算量大而报错，属正常）→ Tasks运行 cg_edge_class_v4_asset。
 * 检查：
 *  - Console“种子点是否落在潮汐候选区”为1；
 *  - 图层“林外：潮汐/鱼塘/陆地”：浅蓝=潮汐区，紫=鱼塘，灰=陆地。前缘泥滩应为浅蓝，基围鱼塘应为紫；
 *    鱼塘被误判为浅蓝 → 画多边形导入为 forceNonTidal；前缘泥滩不是浅蓝 → 增大 CONNECT_MAX_M 或画 forceTidal；
 *  - 图层“林缘类型”：红=向海，橙=向陆，粉=混合，与卫星底图一致。
 * 若导出时报内存不足：把 CONNECT_MAX_M 降到3000再导出。
 ****/

// ======================= 参数 =======================
var aoi = ee.Geometry.Rectangle([106.43, 10.28, 107.36, 10.94]);
var gmwIC = 'projects/sat-io/open-datasets/GMW/extent/GMW_V3';
var PROJ = ee.Projection('EPSG:32648').atScale(10);
var ASSET_ROOT = 'projects/mangrove-494808/assets/'; // GEE资产目录（站点配置 asset_root）

var GMW_SEARCH_PX = 3;        // 在GMW并集外扩3像元(30 m)范围内重新判定红树林
var NDVI_MIN = 0.45;          // 红树林像元：全年NDVI中位数下限（对照底图调整）
var MANGROVE_WATER_MAX = 0.5; // 红树林像元：全年光学水体频率上限
var MIN_PATCH_PX = 9;         // 小于9像元(0.09 ha)的红树林碎斑剔除
var MAX_HOLE_PX = 4;          // 林内非潮汐小空洞(<=4像元)视为冠层缝隙，并入红树林

var WATER_MIN = 0.05;         // 林外像元：全年水体频率>=5%视为可能受潮水影响
var MUD_NDVI_MAX = 0.25;      // 裸潮滩：NDVI中位数上限
var MUD_MNDWI_MIN = -0.25;    // 裸潮滩：MNDWI中位数下限（区分湿泥与干燥建筑/道路）
var MUD_B11_MAX = 1500;       // 裸潮滩：B11中位数上限（湿泥SWIR低，建筑SWIR高）
var BAY_SEED = ee.Geometry.Point([107.07, 10.34]);  // 湾内开阔水面种子点
// 连通起点：种子点 + 自动识别的“开阔湾面”（常年水体且距最近非水体>250 m；鱼塘宽度远小于500 m，不会入选）
var OPEN_WATER_FREQ = 0.9, OPEN_WATER_CORE_M = 250;
var CONNECT_MAX_M = 4000;     // 从开阔湾面向外最多连通的距离（40 km会报 Reprojection output too large）

var EDGE_M = 30, TRANS_M = 60, CLEAR_M = 60;

// 手工修正：在Code Editor左上角画多边形，导入后命名为以下变量，并删掉下面两行同名定义
var forceTidal = ee.FeatureCollection([]);     // 应为潮汐区却未连通的地方（如被云影打断的潮沟）
var forceNonTidal = ee.FeatureCollection([]);  // 误连通的鱼塘/基围
// =====================================================

Map.setCenter(106.895, 10.61, 14);
Map.setOptions('SATELLITE');

// 按波段顺序重命名：直接导出的资产与“Drive导出GeoTIFF再上传”的资产（波段名b1..b7）都能读取；-9999为无数据
var BASE_BANDS = ['ndvi_med', 'mndwi_med', 'b11_med', 'water_freq_all', 'valid_all', 'water_freq_ond', 'valid_ond'];
var base = ee.Image(ASSET_ROOT + 'cg_base_v4').rename(BASE_BANDS);
base = base.updateMask(base.neq(-9999));
var ndviMed = base.select('ndvi_med'), mndwiMed = base.select('mndwi_med'), b11Med = base.select('b11_med');
var waterFreqAll = base.select('water_freq_all');

// ---------- 1. 红树林边界 ----------
var gmw = ee.ImageCollection(gmwIC).max().gt(0).unmask(0).clip(aoi).reproject(PROJ);
var search = gmw.focalMax({radius: GMW_SEARCH_PX, units: 'pixels', kernelType: 'circle'});
// GMW外只允许向“水/泥滩一侧”扩展（30 m内有水体或湿泥）：补上GMW漏掉的向海幼林；
// 向陆一侧的高NDVI像元多为陆生乔灌草，不并入（第一次运行时向陆林缘62%的样点落在GMW外）。
var mudLike = ndviMed.lte(MUD_NDVI_MAX).and(mndwiMed.gte(MUD_MNDWI_MIN)).and(b11Med.lte(MUD_B11_MAX));
var nearWet = waterFreqAll.gte(WATER_MIN).or(mudLike).unmask(0)
  .focalMax({radius: GMW_SEARCH_PX, units: 'pixels', kernelType: 'circle'});
var allowed = gmw.or(search.and(nearWet));
var mang0 = allowed.and(ndviMed.gte(NDVI_MIN)).and(waterFreqAll.lte(MANGROVE_WATER_MAX)).unmask(0);
mang0 = mang0.and(mang0.selfMask().connectedPixelCount(MIN_PATCH_PX + 1, true).gte(MIN_PATCH_PX).unmask(0));
var nonM0 = mang0.not().clip(aoi);

// ---------- 2. 林外潮汐连通区 ----------
var mud = ndviMed.lte(MUD_NDVI_MAX).and(mndwiMed.gte(MUD_MNDWI_MIN)).and(b11Med.lte(MUD_B11_MAX));
var tidalCand = nonM0.and(waterFreqAll.gte(WATER_MIN).or(mud)).unmask(0)
  .where(ee.Image().paint(forceTidal, 1).unmask(0).and(nonM0), 1)
  .where(ee.Image().paint(forceNonTidal, 1).unmask(0), 0);
var notOpenWater = waterFreqAll.lt(OPEN_WATER_FREQ).unmask(1);
var openCore = notOpenWater.fastDistanceTransform(64, 'pixels', 'squared_euclidean').sqrt().multiply(10)
  .unmask(9999).gt(OPEN_WATER_CORE_M).and(notOpenWater.not()).reproject(PROJ);
var seedImg = ee.Image().toByte().paint(BAY_SEED.buffer(50), 1).unmask(0).or(openCore)
  .and(tidalCand).selfMask().reproject(PROJ);
var reach = ee.Image.constant(1).updateMask(tidalCand)
  .cumulativeCost({source: seedImg, maxDistance: CONNECT_MAX_M});
var tidal = reach.mask().unmask(0).and(tidalCand).reproject(PROJ);

// 填林内小空洞（非潮汐、面积很小），避免冠层缝隙造成伪林缘
var smallHole = nonM0.and(tidal.not()).selfMask()
  .connectedPixelCount(MAX_HOLE_PX + 1, false).lte(MAX_HOLE_PX).unmask(0);
var mangrove = mang0.or(smallHole).rename('mangrove').reproject(PROJ);
var nonM = mangrove.not();
var nonTidal = nonM.and(tidal.not());
var outsideClass = ee.Image(0).where(tidal, 1)
  .where(nonTidal.and(waterFreqAll.gte(0.5)), 2)
  .where(nonTidal.and(waterFreqAll.lt(0.5)), 3).updateMask(nonM).rename('outside_class');

// ---------- 3. 欧氏距离（m）----------
function distTo(src) {  // 超出256像元窗口的补9999 m，否则“远离潮汐区”的向陆林缘会被误判为混合
  return src.unmask(0).reproject(PROJ).fastDistanceTransform(256, 'pixels', 'squared_euclidean')
    .sqrt().multiply(10).unmask(9999).reproject(PROJ);
}
var dB = distTo(nonM).rename('d_boundary');
var dT = distTo(tidal).rename('d_tidal');
var dN = distTo(nonTidal).rename('d_nontidal');

// ---------- 4. 林缘类型 ----------
var isEdge = dB.lte(EDGE_M);
var edgeClass = ee.Image(0)
  .where(dB.gt(EDGE_M).and(dB.lte(TRANS_M)), 1)
  .where(isEdge, 4)
  .where(isEdge.and(dT.lte(EDGE_M)).and(dN.gt(CLEAR_M)), 2)
  .where(isEdge.and(dN.lte(EDGE_M)).and(dT.gt(CLEAR_M)), 3)
  .updateMask(mangrove).rename('edge_class');

var g30 = gmw.focalMin({radius: 3, units: 'pixels', kernelType: 'square'}).gt(0);
var g60 = gmw.focalMin({radius: 6, units: 'pixels', kernelType: 'square'}).gt(0);
var gmwZone = ee.Image(0).where(g30.not(), 2).where(g30.and(g60.not()), 1)
  .updateMask(gmw.selfMask()).unmask(-1).rename('gmw_edge_zone');

// ---------- 检查与导出 ----------
print('种子点是否落在潮汐候选区（必须为1）', tidalCand.reduceRegion({reducer: ee.Reducer.max(),
  geometry: BAY_SEED.buffer(50), crs: PROJ}));
Map.addLayer(ndviMed, {min: 0, max: 0.9, palette: ['8c510a', 'f6e8c3', '1a9850']}, 'NDVI全年中位数', false);
Map.addLayer(waterFreqAll, {min: 0, max: 1, palette: ['ffffff', '2166ac']}, '全年水体频率', false);
Map.addLayer(gmw.selfMask(), {palette: ['ffff00']}, 'GMW并集', false);
Map.addLayer(seedImg, {palette: ['ff0000']}, '连通起点（开阔湾面）', false);
Map.addLayer(outsideClass, {min: 1, max: 3, palette: ['00bfff', '8000ff', 'a0a0a0']}, '林外：潮汐/鱼塘/陆地');
Map.addLayer(edgeClass, {min: 0, max: 4, palette: ['1a9641', 'a6d96a', 'd7191c', 'ff8c00', 'ff00ff']},
  '林缘类型 0核心 1过渡 2向海 3向陆 4混合');
Map.addLayer(BAY_SEED, {color: 'red'}, '湾内种子点');

// 资产导出要求各波段类型一致，统一转为float
var out = edgeClass.unmask(255).toFloat().addBands([dB.toFloat(), dT.toFloat(), dN.toFloat(),
  outsideClass.unmask(0).toFloat(), gmwZone.toFloat()]);
// 路线①：直接导出资产
Export.image.toAsset({image: out, description: 'cg_edge_class_v4_asset', assetId: ASSET_ROOT + 'cg_edge_class_v4',
  region: aoi, crs: 'EPSG:32648', scale: 10, maxPixels: 1e10});
// 路线②（资产导出失败时用）：GeoTIFF到Drive → 下载 → Assets面板 NEW → Image upload，资产名填 cg_edge_class_v4
Export.image.toDrive({image: out, description: 'cg_edge_class_v4_drive', fileNamePrefix: 'cg_edge_class_v4',
  region: aoi, crs: 'EPSG:32648', scale: 10, maxPixels: 1e10});
