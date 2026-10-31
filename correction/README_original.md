# 潮位对红树林遥感反演影响：站点流程（交接版）

研究问题：Sentinel-2过境时的潮位，会不会影响红树林的光谱、制图、生物量/碳储量估算和变化检测？影响有多大？
这套流程已在**深圳湾/后海湾**和**广西山口**两个站点完成。现在用于在**潮差更大的地点**重复检验。

---

## 1. 已有结论（新站点的结果可以和这里对照）

| 结论 | 后海湾 | 广西山口 |
|---|---|---|
| FES2022潮时偏差（影像水边线检验） | 提前约2 h | 提前约1.75 h |
| 过境时海水能否到达红树林前缘 | 几乎不能 | 经常能（35/122个日期） |
| 向海林缘MNDWI潮位效应（每米，相对核心区） | +0.033 | +0.040 |
| 向海林缘NDVI潮位效应 | ≈0（不显著） | ≈0（不显著） |
| 向陆林缘（对照） | ≈0 | ≈0 |
| 高潮合成使红树林面积少算（NDVI＋MNDWI去水体制图） | −1.5%（向海林缘−8%） | −6.7%（向海林缘−22%） |
| 生物量密度偏差（固定范围） | ±1.3%以内 | ±1.3%以内 |
| 两期变化检测：后期潮位高≥1 m造成的虚假面积变化 | 约−2.6个百分点 | 约−9个百分点（≈10年真实变化） |

一句话：**NDVI对潮位稳健；对水敏感的波段、水体掩膜制图和林缘像元会产生系统偏差，在变化检测中尤其明显；潮差越大，偏差越大。**
潮差更大的站点，预期效应会更强。如果不是这样，就很值得深究原因。

---

## 2. 目录结构

```
tide_mangrove_pipeline/
├─ README.md                本文件
├─ requirements.txt         Python依赖
├─ pipeline.py              主程序（gee / tides / dates / lagcheck / analyze 五个子命令）
├─ gee_templates/           GEE脚本模板（不要直接运行，由 pipeline.py gee 填入站点参数后生成）
├─ py/                      分析模块（pipeline.py 自动调用）
├─ sites/                   站点配置；每个站点的工作目录也建在这里
│   ├─ example_shankou.json
│   └─ example_houhaiwan.json
├─ tools/hko_tide_match.py  （可选）香港站点：用香港天文台验潮站潮位代替FES
└─ dev/build_templates.py   （开发用）模板由原始后海湾脚本自动转换而来，保留转换规则以便追溯
```

每个站点的工作目录（例如 `sites/mysite/`）由程序自动建立：

```
gee/      生成的GEE脚本（复制到 code.earthengine.google.com 运行）
inputs/   从Google Drive下载的GEE导出文件，保持原文件名，直接放进来
dates/    潮位、合格日期、高/低潮配对日期
results/  分析结果（1_waterline、2_edge_effects、3_tide_composite、4_change_sim）
```

---

## 3. 准备工作（只做一次）

1. **Google Earth Engine账号**，并有一个可写的资产目录（Code Editor左侧Assets面板能看到，例如 `projects/xxx/assets/`）。
   可以先在Assets面板新建一个文件夹，确认资产目录可写。
2. **Python环境**（建议3.11或3.12，单独建一个虚拟环境）：
   ```
   python -m venv .venv
   .venv\Scripts\python -m pip install -r requirements.txt
   ```
3. **FES2022潮汐模型文件**（约5 GB，34个分潮的 `*_fes2022.nc`），从AVISO申请下载，按eo-tides要求的结构放置：
   ```
   <某目录>/fes2022b/ocean_tide_20241025/m2_fes2022.nc ……
   ```
   配置文件里的 `tide.fes_dir` 填 `<某目录>`。
   **该路径不能含中文**（netCDF4在Windows下打不开含中文的路径），或者填相对于运行目录的相对路径。
4. 所有 `python pipeline.py ...` 命令都在 `tide_mangrove_pipeline/` 目录下运行。

---

## 4. 新站点配置

复制 `sites/example_shankou.json`，改名为例如 `sites/mysite.json`，修改以下几项（以下划线开头的键只是说明，可以留着）：

| 参数 | 含义与选法 |
|---|---|
| `site_name`、`prefix` | 站点名称；资产与文件名前缀（如 `zj_`，英文字母加下划线，不要与已有站点重复） |
| `asset_root` | 你的GEE资产目录 |
| `aoi` | 研究区 [西, 南, 东, 北]。包住全部红树林和外侧潮滩，四周留1–2 km。**不宜过大**（建议不超过约30 km×30 km），否则GEE容易超内存 |
| `bay_seed` | 一个落在开阔海面上的点。GEE分类脚本会检查它（Console里必须显示1） |
| `shore_split` | 把红树林分成两侧做分层采样：`axis` 填 `latitude` 或 `longitude`，`value` 为分界线；没有明显两侧时，取研究区中线即可 |
| `tide.lon/lat` | FES潮位代表点，取研究区外侧开阔水面上的点 |
| `ndvi_min` | 红树林判定的NDVI下限（默认0.45）；低矮或稀疏的红树林可调到0.35–0.4 |
| `connect_max_m` | 潮滩与开阔海面的最大连通距离（默认4000 m）。潮滩很宽的地方可增大到6000–8000；GEE报“Reprojection output too large”时减小 |
| `high_tide_threshold_m` | 阈值模型的“高潮日”界限，可在水边线结果出来后参考前缘0–10 m的高程再定 |
| `change_sim` | 两期变化模拟的前期截止年份、后期起始年份 |

---

## 5. 运行步骤

“本地”指在 `tide_mangrove_pipeline/` 目录下运行命令；“GEE”指把 `sites/<站点>/gee/` 里对应的脚本全文复制到GEE Code Editor运行，
然后在右侧 **Tasks** 面板点每个导出任务的 **RUN**。导出的CSV从Google Drive下载后放进 `sites/<站点>/inputs/`。

| 步骤 | 在哪 | 做什么 | 产出 / 检查 |
|---|---|---|---|
| 0 | 本地 | `python pipeline.py gee sites/mysite.json` | 生成9个GEE脚本（日期相关的脚本暂时为空，运行会报“日期尚未写入”，属正常） |
| 1 | GEE | `*_g01_scene_list.js` | 导出 `*_scene_list_fullyear.csv` → inputs/。Console会打印研究区面积、瓦片和红树林面积，核对量级 |
| 2 | GEE | `*_g02_base_layers.js` | 导出资产 `*_base_v4`（10–30分钟）。**等资产出现在Assets面板后再做下一步** |
| 3 | GEE | `*_g03_classify_edges.js` | 先检查：种子点检查值=1；图层“林外：潮汐/鱼塘/陆地”中，前缘滩涂为浅蓝、养殖塘为紫；“林缘类型”与卫星底图吻合。然后导出资产 `*_edge_class_v4` |
| 4 | GEE | `*_g04_sample_points.js` | 导出 `*_edge_sample_points_v4.csv` → inputs/，以及点资产 `*_edge_points_v4` |
| 5 | 本地 | `python pipeline.py tides sites/mysite.json` | 计算每景过境潮位，以及前后±4 h的潮位（几分钟） |
| 6 | 本地 | `python pipeline.py dates sites/mysite.json` | 筛选合格日期，并把日期写进GEE脚本 g06、g07 |
| 7 | GEE | `*_g06_waterline.js`（可同时运行可选的 `*_g05_phenology.js`） | 导出 `*_front_points_daily.csv`、`*_front_bins_daily.csv` → inputs/ |
| 8 | 本地 | `python pipeline.py lagcheck sites/mysite.json` | 用影像水边线检验FES潮时，自动确定时间偏移，生成校正后日期和高/低潮配对日期，并写进 g08、g09 |
| 9 | GEE | `*_g07_spectra.js`、`*_g08_tide_composite_experiment.js`、`*_g09_change_detection_sim.js`（可同时运行） | 导出 `*_edge_spectra_v4_fullyear.csv`、`*_s15_summary.csv`、`*_s15_pixels.csv`、`*_s16_daily_maps.csv` → inputs/ |
| 10 | 本地 | `python pipeline.py analyze sites/mysite.json` | 全部分析，结果在 results/ |

依赖关系：`g02 → g03 → g04`；`tides → dates → g06 → lagcheck → g08、g09`；`g04 + dates → g07`。
`analyze` 会跳过还没有输入文件的分析，所以可以分批运行。

---

## 6. 结果怎么看

| 目录 | 内容 | 重点看 |
|---|---|---|
| `results/fes_lag_check.csv` | FES时间偏移扫描 | 最佳偏移；过境时刻（+0.0h）与最佳偏移的相关差多少。最佳偏移落在±4 h边界时，说明偏差可能更大，需要扩大 `py/s14_fes_lagscan.py` 里的 `LAGS` 范围 |
| `results/1_waterline/` | 水边线与前缘高程 | `WATERLINE_SUMMARY.md`：前缘0–10 m被淹的日期数、前缘高程；`front_wet_vs_tide.png` 应呈S形 |
| `results/2_edge_effects/results.csv` | 光谱面板模型（点＋日期固定效应＋林缘类型×月份） | `main_all` 中 `x_seaward`（向海−核心）与 `x_landward` 的 NDVI、MNDWI、NIR_B8、SWIR_B11；推断看 `p_t_Gmin`、`p_wcr_date`、`p_wcr_twoway` |
| `results/3_tide_composite/S15_RESULTS.md` | 高潮合成 vs 低潮合成：制图面积、生物量、碳储量 | 各制图方法的“全林”和“向海林缘”偏差；固定范围生物量偏差 |
| `results/4_change_sim/CHANGE_SIM_RESULTS.md` | 两期变化检测模拟 | 每米潮位的面积偏差与真实年变化的比值；“后期高≥1 m”组与“潮位相近”组的差值 |

---

## 7. 常见问题

| 现象 | 原因与处理 |
|---|---|
| GEE第一行就报 `Cannot set properties of undefined (setting 'toString')` | GEE网页没加载完整（网络问题），不是脚本错误。Ctrl+F5强制刷新，先运行 `print(ee.Number(1).add(1))` 确认显示2 |
| `XMLHttpRequest` 加载失败 | 网络不稳。换节点后刷新 |
| `User memory limit exceeded` | 交互模式算不动。图层放大到局部看；直接在Tasks里导出（后台限制宽松得多） |
| `Reprojection output too large` | 减小 `connect_max_m` |
| 资产导出失败（对话框显示 `[object Object]`） | 先照常点RUN，看Tasks的状态；仍失败就运行同一脚本里的 `*_drive` 导出任务，下载GeoTIFF后在Assets面板用 NEW → Image upload 上传，资产名与原名相同（下游脚本兼容这种上传的资产） |
| 资产已存在，导出报错 | 重跑前先在Assets面板删除旧资产 |
| `tides` 报找不到 `*.nc` | FES路径含中文或路径错误；按第3节放到纯英文路径 |
| 分类后GMW范围被大量剔除 | 先检查被剔除的像元：若全年NDVI中位数很低（滩涂），属正常；若是稀疏幼林，调低 `ndvi_min` |
| 两期模拟里有离群景 | 程序已按“核心区判为红树林的比例≥0.8”剔除雾霾或薄云影像，并打印被剔除的日期，需人工核对 |
| 香港站点 | 可用 `tools/hko_tide_match.py` 下载香港天文台验潮站逐时天文潮，代替FES（需自行改为站点代码） |

---

## 8. 方法要点与局限

- **林缘分类**：在GMW并集外扩30 m的范围内，用Sentinel-2全年NDVI重新确定红树林边界（只允许向水一侧扩展）。
  林外区分“与外海连通的潮汐区”和“鱼塘/陆地”。用欧氏距离分出核心（>60 m）、过渡（30–60 m）、向海林缘、向陆林缘、混合林缘。
- **潮位**：FES2022在浅湾潮时偏早；用影像水边线做时间偏移扫描来校正（已在后海湾用验潮站数据验证）。
- **推断**：点＋日期双向固定效应；日期×区块双向聚类；wild cluster bootstrap。双向bootstrap偏保守，应报告多种推断结果。
- **生物量参考**：Simard等全球红树林AGB（GEE社区目录），只用来让模型形式接近实际做法，不代表真实生物量；
  Sentinel-2对它的解释度不高（R²约0.1–0.3），所以密度偏差可能被低估。
- 互花米草等非红树植被可能混入向海林缘，需要结合当地资料核查。
- 只评估了遥感观测偏差，不涉及潮汐对真实固碳的影响。

---

## 9. 与原始脚本的对应关系

| 本流程 | 原始脚本（`tide data/` 与 `tide data/code/`） |
|---|---|
| g01 | s01_s02_scene_list_fullyear_gee.js |
| g02 / g03 / g04 | s06a_base_layers_v4 / s06b_classify_v4 / s06c_sample_v4 |
| g05 | s10_phenology_by_edge_class |
| g06 + py/s11 | s11_front_waterline + s11_waterline_elevation.py |
| g07 + py/s13 | s07_extract_edge_spectra_v4 + s13_edge_v4_fullyear.py（依赖 s09、s08） |
| g08 + py/s15 | s15_template + s15_summarize.py |
| g09 + py/s16 | s16_template + s16_change_simulation.py |
| tides / lagcheck | s14_fes_lagscan.py |
| dates | prepare_hhb_dates.py（改名为 prepare_dates.py） |
