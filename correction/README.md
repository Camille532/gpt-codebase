> 本目录为实测潮位修正版。下载方法：返回仓库首页，点击 **Code → Download ZIP**，解压后使用 `correction/` 内的文件。原项目代码位于仓库根目录。

# 将 FES2022 替换为 UHSLC 实测潮位：修改与最小重跑说明

本版基于本次上传的独立脚本修改。原流程留在 `README_original.md`。原说明中的 `pipeline.py`、`sites/`、模板生成器没有随附件提供，本次不依赖这些缺失文件，下面命令直接运行现有脚本。

## 先明确本次潮位文件

`h383 (1).csv` 无表头，五列依次为 `year, month, day, hour, sea_level_mm`。用户已确认时间为 UTC。UHSLC 的 h383 元数据确认单位为毫米、`-32767` 为缺测；383 对应越南头顿 Vung Tau 站，官方站点表基准为 `station zero`。

- [h383 官方数据说明：单位、缺测值、质量级别](https://uhslc.soest.hawaii.edu/opendap/fast/hourly/h383.nc.info)
- [UHSLC 逐小时站点表：Vung Tau、383、station zero](https://uhslc.soest.hawaii.edu/erddap/tabledap/global_hourly_fast.subset)
- [原下载入口](https://uhslc.soest.hawaii.edu/data/)

程序保留站点零点，只除以1000换算成米，不自动减均值，不把站点零点称作平均海平面。文件共165,492小时，2007-10-15 12:00 UTC至2026-08-31 23:00 UTC，其中4,401小时缺测；2020年8,784小时中只有5,900小时有效。逐年覆盖统计和原g07分析日期的缺测检查见 `TIDE_FILE_AUDIT.md`。

h383 是潮位站实测水位，和 FES 天文潮的暴露定义有差别，可能包含非天文海面变化；它代表头顿站水位，不自动等于研究区每一个林缘点的当地水位。保留水边线一致性诊断，但不把旧 FES 的最佳时间偏移施加到实测数据，也不根据 NDVI 相关性自动挑选偏移。

## 按你提供的芹椰目录直接运行

已按 `D:/桌面/Vietnam/sites/cangio/inputs/cg_scene_list_fullyear.csv` 设置入口 `run_cangio_observed.py`。把修改后的Python文件放在同一目录运行；原项目目录和原输出均保留。

1. 将 `h383 (1).csv` 放进 `D:/桌面/Vietnam/sites/cangio/inputs/`。若已在别处，可加 `--station-file "实际路径"`，无需复制。
2. 查看脚本发现的原日期表（仅读取，不计算）：

```bash
python run_cangio_observed.py --list-dates
```

3. 运行过境匹配、更新原日期表副本、更新GEE分组，并对存在原导出CSV的步骤运行本地分析：

```bash
python run_cangio_observed.py --run-local
```

所有新文件在 `D:/桌面/Vietnam/sites/cangio/observed/`。入口扫描原 `dates/`，优先匹配g07日期集合及原流程的corr最终分析表，记录实际选中的路径；有歧义会列出候选，不会根据文件修改时间猜选。通过 `--analysis-dates "实际文件路径"` 可显式指定；`--dates-dir` 和 `--gee-dir` 可指定不同目录。

若原日期在实测文件中缺测，程序先停止并列出日期；只有接受缩小样本后才加 `--drop-unmatched`：

```bash
python run_cangio_observed.py --drop-unmatched --run-local
```

入口仅为原导出文件确实存在的s11、s13和s16 R1/R2生成/执行命令，具体见 `observed/migration_inputs.json`。缺失导出文件不会被伪造。新GEE脚本在 `observed/gee/`，按其中的 `rerun_plan.json` 决定g08/g09重跑；入口不执行GEE导出，也不自动把旧RF结果作为新低潮训练结果。

以上入口沿用s13默认bootstrap 1999次及s16原默认年份/质量设置。若你原先改过参数，使用后文分步命令明确传入原参数。

## 哪些地方改了、哪些结果可复用

| 文件/步骤 | 本次修改 | 已有结果如何处理 |
|---|---|---|
| `s14_fes_lagscan.py` | 新增 `--station-file`、UHSLC无表头预设、UTC匹配、毫米转米、缺测与插值审计；FES模式保留 | 必须重建场景过境潮位；实测模式无需eo-tides或FES文件 |
| `tide_data.py`（新增） | 统一读表、单位/基准记录、时间匹配、日潮位校验 | 各分析脚本共用，须放在同一目录 |
| `refresh_observed_dates.py`（新增，推荐） | 更新旧日期表的潮位与来源，保留原日期顺序、质量统计、`qualified_strict`等筛选标志；清除旧FES lag列 | 可分别更新“全部质量日期表”和“实际分析日期表”，无需重新做影像质控 |
| `prepare_dates.py` | 可读取任意来源过境潮位；修复 `--required-tiles` 原来被硬编码4覆盖；新建日期表输出严格筛选标志 | 仅从头建日期表时运行；已有结果优先使用refresh，避免更改原分析样本 |
| `cg_g01`—`cg_g05` | 仅g01潮位说明更新，光学计算不变 | 资产、分类、采样点、物候导出可复用；g01清单可直接复用 |
| `cg_g06_waterline.js` | 修正基准说明 | 目标日期已经导出则直接复用；缺少目标日期才补导出 |
| `cg_g07_spectra.js` | 光学提取算法未改；日期由同步工具更新 | 目标日期已经导出则复用，光谱数值无需因换潮位而重算 |
| `s08` / `s09` / `s13` | 从更新日期表读取潮位并校验日期唯一、潮位完整；s13取消默认0.8m阈值 | 重跑你原来实际使用的模型。全年v4分析通常重跑s13即可，不必额外重复s08/s09 |
| `s11_waterline_elevation.py` | 高程、图轴和报告沿用输入潮位基准；动态日期数；修正重复潮位不可拆分的阈值候选 | 原水边线影像表可复用，但高程阈值、湿润比例与潮位关系必须本地重算 |
| `make_tide_composite_dates.py` | 用新潮位重排同月高低潮组；保留基准说明，稳定排序，防止相同潮位虚造对比 | 重新生成分组 |
| `sync_gee_tide_dates.py`（新增） | 自动替换g06/g07/g08/g09日期数组，生成独立脚本副本与新旧分组比较 | 只改日期；AOI、资产、分类、采样和光学公式保留 |
| `cg_g08` + `s15` | g08内置旧FES分组有运行保护；同步工具生成实测分组后解除保护 | 高或低组成员改变就重跑g08，再重跑s15；分组完全相同可复用原GEE导出 |
| `cg_g09` + `s16` | g09日期由同步工具更新；s16新增 `--methods R1 R2` | R1/R2逐日统计不依赖潮位，可先复用并本地重算；低潮训练组改变时RF必须重新训练并导出所有目标日期 |

`s15_summarize.py` 本身不读取潮位，计算逻辑无需修改。它消费g08的分组结果，不能只改报告里的组名就当作实测分组结果。

## 推荐：保留原样本，最少重跑

所有命令从脚本所在目录执行。以下 `old/...`、`inputs/...` 是占位路径，请换成你已有文件的真实路径。每条命令为一行，Windows PowerShell/cmd也可直接使用；含空格的文件名加双引号。所有实测结果写入新的 `observed/`，不要覆盖旧FES结果。

### 1. 只重新匹配过境潮位

直接读取原g01清单，或含 `overpass_utc` 的旧逐景潮位表；必须保留精确UTC过境时刻，不能用一个固定卫星过境时间代替，也不能用日平均潮位。

```bash
python s14_fes_lagscan.py --scene-list inputs/cg_scene_list_fullyear.csv --station-file "h383 (1).csv" --station-format uhslc --station-timezone UTC --station-name "UHSLC383 Vung Tau" --out observed/cg_scene_tides_observed.csv --allow-unmatched
```

UHSLC预设自动使用毫米、缺测码-32767、站点零点。默认仅在相邻有效记录间隔不超过90分钟时线性插值：逐小时数据中间少一个小时就不会跨过去。超出记录范围不外推，不用FES回填；输出 `tide_match_status`、前后观测时间和配对间隔，并写 `.report.json`。`--allow-unmatched` 只允许输出NaN供审计，不代表允许下游直接分析NaN。

同一天多瓦片沿用原方法取场景潮位均值；只要其中一景缺测，该日就标为缺测，防止用部分瓦片掩盖缺测。若同日过境跨度超过30分钟，日期更新会停止，需按轨道/时刻拆分，而不是平均多个不同时相。

### 2. 更新已有的日期表

将你原来实际使用的全部质量日期表、最终分析日期表都传入。输出保留输入文件名，放在新目录；旧表中的 `month`、`valid_mean`、`qualified_strict`等保留，旧 `tide_lag_*`等被移除。

```bash
python refresh_observed_dates.py --tides observed/cg_scene_tides_observed.csv --existing-dates old/cg_dates_all_quality.csv old/cg_analysis_dates_relaxed06_corr.csv --out-dir observed/dates
```

默认某个原日期缺实测潮位就停止，明确列出日期。检查后若接受仅分析实测覆盖日期，再加 `--drop-unmatched` 运行：

```bash
python refresh_observed_dates.py --tides observed/cg_scene_tides_observed.csv --existing-dates old/cg_dates_all_quality.csv old/cg_analysis_dates_relaxed06_corr.csv --out-dir observed/dates --drop-unmatched
```

输出 `refresh_report.json` 记录删掉的日期，每张表还有 `_tide_comparison.csv`，有相应协变量时也重算潮位混叠Spearman诊断。差值同时含模型差异和基准差异，不应直接当作FES精度误差。若有剔除日期，比较FES与实测模型时也应把FES结果限制到相同日期，区分潮位更换与样本变化的影响。

### 3. 先重跑本地水边线与主模型

```bash
python s11_waterline_elevation.py --points inputs/cg_front_points_daily.csv --bins inputs/cg_front_bins_daily.csv --dates observed/dates/cg_dates_all_quality.csv --analysis-dates observed/dates/cg_analysis_dates_relaxed06_corr.csv --out-dir observed/results/waterline
python s13_edge_v4_fullyear.py --spectra inputs/cg_edge_spectra_v4_fullyear.csv --dates observed/dates/cg_analysis_dates_relaxed06_corr.csv --out-dir observed/results/edge_effects
```

s13默认运行连续潮位模型（含严格样本、GMW子样本、向陆基准），不运行绝对高潮阈值模型。确认阈值与站点零点一致后，才补充 `--high <同基准阈值米值>`；不能直接沿用旧 `0.8`。若你的旧日期表缺 `qualified_strict`，需要从原筛选记录补回该列，不能默认把所有日期当作严格合格。

仅需要连续线性潮位效应和同月高低潮排序时，可以保留站点零点；稳定的常数零点偏移不会改变排序、两期潮位差或含固定效应线性模型的斜率。绝对高程、固定高度阈值、跨站点绝对高度比较则必须统一基准。s08还包含平方项，不能把平方模型中的一次项系数简单当作不受零点平移影响的线性斜率。

如你原来使用s08或s09，用同一新日期表另行重跑：

```bash
python s08_analyze.py --spectra inputs/your_existing_spectra.csv --dates observed/dates/cg_analysis_dates_relaxed06_corr.csv --out-dir observed/results/s08
python s09_edge_datefe_bootstrap.py --spectra inputs/your_existing_spectra.csv --dates observed/dates/cg_analysis_dates_relaxed06_corr.csv --out-dir observed/results/s09
```

### 4. 更新GEE日期，并看重跑清单

```bash
python make_tide_composite_dates.py --dates observed/dates/cg_analysis_dates_relaxed06_corr.csv --site "越南芹椰" --out-js observed/high_low_dates.js
python sync_gee_tide_dates.py --dates observed/dates/cg_analysis_dates_relaxed06_corr.csv --waterline-dates observed/dates/cg_dates_all_quality.csv --gee-dir . --out-dir observed/gee
```

生成 `observed/gee/cg_g06...` 至g09四份脚本，以及 `rerun_plan.json`。查看g08分组是否改变、g09的LOW_DATES是否改变。g06/g07数组删除日期时不必重导已有影像表，日期表的本地合并会筛掉它们；若新增日期才补导出。

发布包g08/g09中的原日期只是旧FES日期，不是已计算好的实测日期，所以在同步前会主动报错。运行同步后使用 `observed/gee/` 的脚本。

g08分组改变：到GEE运行新g08，下载新 `*s15_summary*.csv`，放到单独目录后运行：

```bash
python s15_summarize.py --dir observed/exports/s15
```

不要把新旧s15汇总CSV放进同一个目录，否则glob会一起读入。

### 5. 两期变化：先复用R1/R2，RF按训练组决定

```bash
python s16_change_simulation.py --sites "越南芹椰|inputs/cg_s16_daily_maps.csv|observed/dates/cg_analysis_dates_relaxed06_corr.csv" --methods R1 R2 --out-dir observed/results/change_rules
```

`--early-end`、`--late-start`、`--min-cover`、`--qc-core-r1` 保持你原来用的值；示例未写时仍沿用脚本原默认值。

如果LOW_DATES改变，在GEE重跑新g09（包含重新训练RF和应用到全部目标日期），下载新逐日统计，再跑三种方法：

```bash
python s16_change_simulation.py --sites "越南芹椰|observed/exports/cg_s16_daily_maps.csv|observed/dates/cg_analysis_dates_relaxed06_corr.csv" --methods R1 R2 RF --out-dir observed/results/change_all
```

如果LOW_DATES相同、目标日期也已导出，可以直接用旧逐日统计跑三种方法。本地的潮位回归、潮位差分组、配对统计仍必须重算。

## 如果确实要从头重建日期表

```bash
python prepare_dates.py --scene-list inputs/cg_scene_list_fullyear.csv --tides observed/cg_scene_tides_observed.csv --min-valid 0.6 --strict-min-valid 0.8 --required-tiles 4 --out-dir observed/rebuilt_dates
```

这里只是示例，不表示你研究区必须4瓦片或必须采用0.6/0.8；按原站点影像覆盖和既定质量标准填写。重建方式可能改变样本，和“只替换潮位”的结果要区别保存。输出的 `dates_all_quality.csv` 可能包含未匹配的非合格日期；用于水边线分析前应明确处理这些缺测，不能让NaN进入高程阈值。

## 输入格式与验证范围

若使用自行整理的带表头标准化文件（时间列 `time` 为UTC，潮位列 `tide_height_m` 为米，缺测为空值），其基准仍为站点零点。作为输入时请使用 `--station-format columns --station-time-col time --station-unit m --datum "UHSLC station zero"`，不要再使用uhslc无表头预设。

已用真实h383文件验证读取/转换/插值，用模拟数据验证日期更新、分组同步、s11/s13/s16与阈值模型；Python编译、GEE JavaScript语法检查通过。GEE服务端导出未在本会话执行。没有上传原g01 CSV、原日期表及GEE导出结果，因此尚不能生成真实的新过境潮位表、新高低潮分组或新的研究结论。所有示例路径需换成你的现有文件名，不能拿模拟验证输出当作研究结果。
