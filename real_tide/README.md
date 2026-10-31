# real_tide：按原 Vietnam 工程结构替换 FES2022

本目录基于本次上传的 `Vietnam.zip` 修改。只替换潮位数据来源、撤掉FES时差扫描，并处理更换高程基准后不能照搬的绝对阈值。没有新增实测与模拟潮位比较，没有改动光学指数、分类算法、回归公式、固定效应、bootstrap或原质量筛选规则。

使用原 `pipeline.py` 入口。不要再运行前一版 `run_cangio_observed.py`；前一版根目录的同名独立Python文件不参与本版运行，本版从 `py/` 加载模块。

## 1. 只覆盖这五个文件

GitHub 的 `real_tide/` 只是交付目录，不是电脑上需要新增的运行目录。把其内容按相对路径覆盖到 `D:/桌面/Vietnam/`：

| GitHub文件 | 本机保存位置 | 核心变化 |
|---|---|---|
| `pipeline.py` | `D:/桌面/Vietnam/pipeline.py` | tides读取文件；dates保留原样本质量标志并替换潮位；不再调用FES lag扫描；仍生成原文件名和GEE脚本；阻止用旧分组导出冒充新分组结果 |
| `py/s14_fes_lagscan.py` | `D:/桌面/Vietnam/py/s14_fes_lagscan.py` | 彻底删除eo-tides、FES模型与±4小时计算，只读UHSLC无表头CSV。旧文件名和run_model函数名保留，仅为兼容原入口 |
| `py/s13_edge_v4_fullyear.py` | `D:/桌面/Vietnam/py/s13_edge_v4_fullyear.py` | 连续潮位模型完全不变；绝对高潮阈值改为显式可选，未指定就只跳过这一项 |
| `py/s11_waterline_elevation.py` | `D:/桌面/Vietnam/py/s11_waterline_elevation.py` | 图轴、报告从日期表读取实际高程基准和日期数；水边线高程算法不变 |
| `sites/example_cangio.json` | `D:/桌面/Vietnam/sites/example_cangio.json` | tide.source改为uhslc_file，file指向实测文件，删除FES路径/代表点/时差；高潮绝对阈值设为null |

本目录没有打包原工程的其他脚本；它们直接沿用你原来的版本。安装时先保留这五个原文件的副本。`gee_templates/`、其他 `py/` 文件、已有GEE资产、`inputs/`导出不需要整体替换或搬动。

目录应为：

```text
D:/桌面/Vietnam/
  pipeline.py
  py/
    s14_fes_lagscan.py
    s13_edge_v4_fullyear.py
    s11_waterline_elevation.py
    ...原来的其他Python文件
  gee_templates/                   原模板保留
  sites/
    example_cangio.json
    cangio/
      inputs/
        h383 (1).csv               放入实测文件
        cg_scene_list_fullyear.csv
        ...原来的GEE导出文件
      dates/                       继续保存中间表和GEE日期列表
      gee/                         继续保存生成后的GEE脚本
      results/                     继续保存分析结果
```

站点AOI、资产地址、前缀cg_、2019—2024年份、0.6/0.8质量阈值、bootstrap_B=1999及变化模拟年份均保留上传配置。唯一额外停用的是原相对平均海平面的1.0m绝对阈值，原因见下文。

## 2. 潮位到底怎样替换

文件为你提供的 `h383 (1).csv`，无表头五列：年、月、日、UTC小时、潮位毫米值。

1. 根据年/月/日/小时构建UTC时间，和场景清单的 `overpass_utc` 对齐，不改成北京时间。
2. `-32767`及空潮位记录作为缺测；毫米除1000得到米。
3. 精确命中有效小时就直接用该观测；其余仅在前后相邻有效记录间线性插值。
4. 最大允许间隔90分钟。逐小时数据少一个小时，前后有效记录间隔会变成2小时，这种缺口不跨越；范围外不外推，不用FES填补。
5. 不减均值、不自动换零点、不按水边线或NDVI相关性挑时差。输出基准为 `UHSLC station zero`。
6. 同一天多瓦片仍沿用原流程取场景潮位均值；任一场景缺潮位就将该日标缺测，不用部分瓦片的均值掩盖缺测。

`dates/scene_tides.csv`仍保留原场景字段和 `tide_height_m`，增加匹配状态、前后观测时刻、来源与基准。原FES的 `tide_lag_*`、`tide_fes0`、`tide_corr`和FES时差字段不再传入新的活动日期表。

原来的 `analysis_dates_corr.csv`、`all_dates_corr.csv`和 `analysis_dates_corr_strict.csv` 文件名继续保留，是为了让现有下游脚本直接工作。现在 `_corr` 文件只是实测日期表接口，不代表仍有FES校正。

原 `lagcheck` 命令仅保留为 `dates` 的兼容别名，不做相关扫描、不产生新的fes_lag_check.csv，也不需要运行。旧 `--lag` 参数会被拒绝。

UHSLC依据：[h383单位与缺测值元数据](https://uhslc.soest.hawaii.edu/opendap/fast/hourly/h383.nc.info)、[Vung Tau站383及station zero](https://uhslc.soest.hawaii.edu/erddap/tabledap/global_hourly_fast.subset)。实测站水位代表Vung Tau站，不自动等同所有林缘位置的当地水位；本次按你的要求只替换输入，不拟合当地传播时差或比较FES。

## 3. 已经跑过原流程，现在只运行这些命令

在CMD执行，每条命令一行。使用你原来安装numpy/pandas/matplotlib的Python环境。实测读取不再需要eo-tides或FES文件。

先进入原工程根目录：

```cmd
cd /d "D:\桌面\Vietnam"
```

### 第一步：重新匹配实测过境潮位

```cmd
python pipeline.py tides sites/example_cangio.json
```

输入：`inputs/cg_scene_list_fullyear.csv`、`inputs/h383 (1).csv`。

输出：`dates/scene_tides.csv`、`dates/scene_tides.report.json`。这一步完成只表示匹配完成，模型还没重跑。

### 第二步：更新原日期表及GEE日期分组

你这份数据存在缺测。如果接受仅使用实测完整覆盖日期，运行：

```cmd
python pipeline.py dates sites/example_cangio.json --drop-unmatched
```

如果先想检查缺测，可以不加 `--drop-unmatched`，程序会列出缺测数量后停止，写出 `dates/real_tide_missing_dates.csv`，不会把NaN送进模型。

此处从原 `analysis_dates.csv`和 `prep/dates_all_quality.csv`保留既有日期质量统计与严格标志，只替换潮位并按明确的缺测规则删日；不重新调整瓦片数、有效比例或季节筛选。首次运行会自动在 `dates/_before_real_tide/`保存必要原日期表和旧GEE分组，反复运行仍以原质量样本为起点。

继续输出原文件名：

- `analysis_dates.csv`、`analysis_dates_corr.csv`：新的完整实测分析日期。
- `analysis_dates_corr_strict.csv`：保留原严格标志的实测子样本。
- `all_dates_corr.csv`：所有能匹配实测潮位的质量日期。
- `prep/dates_all_quality.csv`：仍保留全部855个质量日期，缺潮位用NaN记录，仅作为质控中间表。
- `tide_composite_dates.js`：仍使用原make_tide_composite_dates.py的同月高低潮分组算法。
- `sites/cangio/gee/`：按原模板生成更新日期后的脚本；这里只写JS文件，不执行GEE。
- `results/real_tide_status.json`：日期计数、需补导出的日期以及g08/g09是否必须重跑。

这一步完成仍只是潮位/日期/脚本准备完成。

### 第三步：先重跑水边线和光谱模型

```cmd
python pipeline.py analyze sites/example_cangio.json --steps waterline,edge
```

复用 `inputs/cg_front_points_daily.csv`、`cg_front_bins_daily.csv`、`cg_edge_spectra_v4_fullyear.csv`，更新 `results/1_waterline/`、`results/2_edge_effects/`。

这里保留原min-obs=40、原回归公式和bootstrap次数；模型可能运行较久。s08/s09仅作为原s13依赖，无需单独跑。分析输出会写回原results子目录，如需保留旧分析结果，运行前先复制这些旧结果子目录存档；代码不会自动比较旧结果。

### 第四步：在GEE只重跑g08和g09

本次实际数据已确认高低潮分组和低潮训练组都改变，故需运行：

```text
D:/桌面/Vietnam/sites/cangio/gee/cg_g08_tide_composite_experiment.js
D:/桌面/Vietnam/sites/cangio/gee/cg_g09_change_detection_sim.js
```

复制更新后的脚本到GEE Code Editor运行，并在Tasks提交导出。g08/g09可以分别提交。

- g08下载新的 `cg_s15_summary.csv`和 `cg_s15_pixels.csv`，按原名替换 `inputs/`中的同名导出。
- g09下载新的 `cg_s16_daily_maps.csv`，按原名替换 `inputs/`中的同名导出。RF训练组变了，需对全部新目标日期重新应用RF；只改本地潮位列不能更新这个RF模型。

g08/g09的分类、训练参数、光谱、生物量公式均保持原模板；变化仅为重新填入实测筛出的日期。

### 第五步：新GEE导出下载完，再做最后两项本地分析

```cmd
python pipeline.py analyze sites/example_cangio.json --steps composite,change
```

分别更新 `results/3_tide_composite/`与 `results/4_change_sim/`。如果inputs里的文件仍与原导出内容相同而GEE分组已变，入口会明确跳过相应步骤，防止误把旧分组统计当作新实测结果；这种“跳过”不代表整个任务已完成。

如果重新导出的文件碰巧内容逐字节完全相同，也会触发上述保护，此时请先核实新导出任务确实使用了本次新脚本，再反馈，不要更改CSV内容来绕过检查。

## 4. 哪些步骤不用跑，哪些必须跑

| 原步骤 | 本次是否重跑 | 原因 |
|---|---|---|
| g01场景清单 | 不用 | 已有精确UTC过境时间 |
| g02基础图层、g03林缘分类、g04样点 | 不用 | 不读取FES数值，资产和样点沿用 |
| g05物候 | 不用 | 指标计算未改变 |
| g06水边线导出 | 不用 | 已有导出覆盖新99个分析日期 |
| g07光谱导出 | 不用 | 已有导出覆盖新99个分析日期，指数数值不因潮位来源变化而重算 |
| pipeline tides | 必须 | FES计算改为文件匹配 |
| pipeline dates | 必须 | 替换潮位、排除缺测、重新生成高低潮分组 |
| pipeline lagcheck | 不用 | FES扫描和自动时差校正已移除 |
| s11水边线本地分析 | 必须 | 潮位、高程零点和拟合输入改变，干湿观测复用 |
| s13光谱面板模型 | 必须 | 更换模型中的潮位自变量，光谱观测复用 |
| g08高低潮合成 | 必须 | 本次高低潮日期组改变 |
| s15本地汇总 | g08新导出后必须 | 消费新的分组统计 |
| g09变化模拟导出 | 必须 | 本次RF低潮训练组改变 |
| s16本地变化分析 | g09新导出后必须 | 用新的过境潮位和一致的RF统计重新分析 |

本版为保持其他原脚本不变，s16仍按原流程一次分析R1/R2/RF，没有再加入上一版单独先跑R1/R2的分支。

## 5. 高程基准与绝对高潮阈值

UHSLC值除1000后是相对站点零点的米值，不是相对平均海平面。保留这个基准即可做连续潮位模型、同月高低潮排序及两期潮位差分析；不能把原配置相对平均海平面的 `high_tide_threshold_m=1.0`直接套上。

配置现为 `"high_tide_threshold_m": null`：仅跳过 `high_tide_threshold` 模型，其他四组连续潮位模型及其推断保持原样。以后给出明确的、相对同一站点零点的阈值，可填入数值恢复原阈值模型。没有擅自用样本均值或分位数替代旧阈值，也没有拟合FES零点偏移。

s11仅改运行输出中的基准标签，没有改变waterline_threshold计算。模板与部分原脚本头部仍可能有旧FES注释，本次按你的要求没有为纯注释变化打包额外文件；活动数据以tide_source/tide_datum和本说明为准。

## 6. 已用你这份工程验证的事实

用上传的场景清单、原质量表及h383实际匹配：

- 实测原始165,492小时，4,401小时缺测；场景匹配有127景缺测。
- 原全部质量日期855天，46天没有完整实测过境潮位，保留809天供水边线分析接口使用。
- 原分析105天，缺6天，保留99天；原严格48天，保留43天。
- 缺少实测的6个分析日期：2020-11-18、2020-12-03、2021-03-13、2024-01-12、2024-02-11、2024-02-26。
- 原水边线和光谱导出已覆盖新目标日期，g06/g07无需补导出。
- g08分组及g09低潮训练组均改变，需要GEE重跑。

匹配、缺测/插值、日期更新与原函数未变检查通过；用你的真实导出运行了s11/s13的端到端检查。检查时仅把临时测试配置的bootstrap_B设为0以验证流程，交付配置仍为1999，正式结果请按上述命令在本机重跑。未执行实测与FES关系比较，也未在GEE提交导出。这里没有上传潮位CSV、inputs数据或研究结果。

**完成判断：先完成第三步的本地模型，再完成第四步的GEE导出，最后完成第五步的本地汇总，才是本轮全部完成。`type`仅显示报告，不能用来启动或判断模型运行。**
