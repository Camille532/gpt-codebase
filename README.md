# 潮位对红树林遥感反演影响

这份仓库整理了 Sentinel-2、GEE 和 Python 分析脚本，用于研究过境潮位对红树林光谱、林缘分类、制图面积、生物量估算及变化检测的影响。原交接流程保存在 [README_handover.md](README_handover.md)；其中描述的 `pipeline.py`、站点配置和模板目录未随本次代码附件提供，不能直接照其自动化命令运行。

## 当前脚本工作流

1. 在 GEE 运行 `cg_g01_scene_list.js`，导出影像清单 CSV（需包含 `overpass_utc`）。
2. 将潮位站数据准备为 CSV，含时间列和潮位（米）列。运行 `s14_fes_lagscan.py --station-csv ...`，按 UTC 过境时间插值并生成带 `tide_height_m` 的场景潮位表。站点无时区标记时，必须用 `--station-timezone` 指定其真实时区；程序将其转换为 UTC。时间不能覆盖过境时刻或观测间隔超过 `--max-gap` 时程序会报错，不做外推。
3. 运行 `prepare_dates.py --scene-list ... --tides ... --out-dir ...`，生成日期质量表。`--required-tiles` 和 `--min-valid` 控制筛选条件。
4. 在 GEE 按依赖运行 `cg_g02_base_layers.js` → `cg_g03_classify_edges.js` → `cg_g04_sample_points.js`。各脚本输出的资产需成功创建后再运行下游脚本。
5. 运行 GEE `cg_g05_phenology.js`（可选）及 `cg_g06_waterline.js`，导出物候、水边线数据；本地运行 `s11_waterline_elevation.py` 分析水边线与前缘高程。
6. 使用合格日期表更新/填写 `cg_g07_spectra.js` 的日期后运行，导出光谱表。本地运行 `s08_analyze.py`、`s09_edge_datefe_bootstrap.py` 或 `s13_edge_v4_fullyear.py` 做面板分析。`s13` 需要日期表提供 `qualified_strict`。
7. 运行 `make_tide_composite_dates.py` 生成同月高低潮日期，再将日期填入 `cg_g08_tide_composite_experiment.js` 并运行；用 `s15_summarize.py` 汇总结果。
8. 运行 `cg_g09_change_detection_sim.js` 导出逐日分类统计，再运行 `s16_change_simulation.py` 评估两期变化偏差。

潮位站数据替代 FES 后，FES 时间偏移扫描不再适用。水边线分析仍可独立用于核对影像潮滩湿润情况。各本地脚本的输入输出路径以命令行帮助及源码中的参数为准；大部分分析需要相应 GEE 导出数据。

## 依赖

Python 3.10+，以及 `numpy`、`pandas`；图表输出还需要 `matplotlib`。只有使用 FES2022 模式时才需要 `eo-tides` 与 FES 数据文件。GEE JavaScript 脚本需在 Earth Engine Code Editor 执行。

```bash
python -m pip install numpy pandas matplotlib
```

## 附件整理

仓库包含随任务提供的 9 个 Python 与 9 个 GEE JavaScript 独立脚本。重复上传且内容相同的 GEE 文件已去重。
