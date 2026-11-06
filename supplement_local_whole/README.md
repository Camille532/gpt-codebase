# 红树林潮位研究：局部响应—全林估计补充分析

这套代码用于已有的 Vietnam 项目，也可以直接读取本次 cangio.zip 解压后的目录。原来的 pipeline.py、py/、输入CSV和GEE资产不需要替换。程序把补充结果写到新的结果文件夹。

建议先用 Python 完成统计和汇总，再用 Origin 作图。现有数据就能完成第一轮补充；只有细分过渡带的全区面积统计和共同有效范围碳量分解，需要额外运行 GEE。

## 一、先明确分区，避免研究对象混淆

原 g03 的定义是：

| 原编码 | 定义 | 研究中的用途 |
|---|---|---|
| 0 | 核心区：距最近非红树林像元 >60 m | 光谱模型参照，不代表已证明不受潮位影响 |
| 1 | 过渡带：距最近非红树林像元 >30且≤60 m | 包含不同方向的过渡带 |
| 2 | 向海林缘：边界内30 m，邻近潮汐连通区，并远离非潮汐区 | 主要研究对象 |
| 3、4 | 向陆、混合林缘 | 保留在模型控制和全林总量中 |

**原过渡带不是专门的向海过渡带，向海林缘＋过渡带也不是已确定的潮间带。**

本包额外检查一个几何代理分区：

- 靠潮汐区一侧的过渡带：原编码1，且 `d_tidal < d_nontidal`。
- 其余过渡带：原编码1，且 `d_tidal >= d_nontidal`；距离相等归到这一组。

这里的“靠潮汐区一侧”只表示更靠近已有潮汐连通区。它不证明林下淹没、不模拟水流，也不保证样点处在同一横断面上。新分区基于固定空间信息，不按显著性选点。它属于本次探索性补充，后续应在其他时期/站点验证。

本次数据中原过渡带有300个样点，其中69个符合上述靠潮汐区条件。代码保留完整原分区，并另做拆分后的模型，没有删去其余样点。

## 二、先跑现有数据：无需重跑GEE

### 1. 放置代码

把整个 `supplement_local_whole` 文件夹复制到项目根目录：

```text
D:/桌面/Vietnam/
  pipeline.py
  py/
  supplement_local_whole/
    local_whole.py
    make_supplement_gee.py
    requirements.txt
    test_supplement.py
    README.md
  sites/
    cangio/
      inputs/
      results/
      gee/
```

如果下载的是GitHub整个仓库的ZIP，只取其中这个文件夹放入你的项目；不要把仓库当成一套新的Vietnam项目覆盖。

### 2. 安装依赖并运行

Windows CMD 中运行：

```bat
cd /d "D:\桌面\Vietnam"
python -m pip install -r supplement_local_whole/requirements.txt
python supplement_local_whole/local_whole.py --site-dir sites/cangio
```

完成后结果在：

```text
sites/cangio/results/5_local_whole/
```

程序默认从已有 `results/1_waterline/date_front_wetness.csv` 读取日期潮位。本次归档中该文件的 `analysis_date` 对应严格日期集合，代码据此复现48天的严格模型；主模型按光谱文件实际有数据、且潮位有效的日期取交集。

如果日后更换数据版本，可以明确指定原日期表：

```bat
python supplement_local_whole/local_whole.py --site-dir sites/cangio --dates sites/cangio/dates/analysis_dates_corr.csv
```

日期表需含 `date`、`tide_height_m`，严格模型另需 `qualified_strict`。不要给实测影像分组配FES合成导出，再把它说成同一版本完整实验。

本包读取现有FES潮位，没有重新估算FES，也没有把旧移时结果自动替换成0小时结果。若要比较不同潮位来源，必须另建输出目录，并统一日期、分组及对应GEE导出。

### 3. 五步对应哪些结果

| 研究步骤 | 现有数据可完成什么 | 输出 |
|---|---|---|
| ① 固定分区、核实样本 | 分区参考面积、样点数、日期数 | `00_zone_area_profile.csv`、`00_date_audit.csv`、`01_zone_samples.csv` |
| ② 局部光谱响应 | 向海林缘、过渡带相对核心区的斜率；两者直接差异；拆分过渡带后的补充模型 | `01_spectral_contrasts.csv` |
| ③ 局部制图与生物量 | 分区高低潮面积、固定范围生物量、随制图范围汇总的地上生物量及碳量 | `02_composite_local_whole.csv` |
| ④ 局部对全林的贡献 | 向海林缘＋过渡带、其余森林、全林的绝对差异及共同分母下的贡献 | `02_composite_local_whole.csv`、`03_contribution_checks.csv` |
| ⑤ 跨期监测应用 | 各分区逐日面积模型、局部覆盖率敏感性分析、前后期配对 | `04_daily_area_and_qc.csv`、`05_daily_area_effects.csv`、`06_pair_summary.csv` |

`READ_RESULTS.md` 说明统计口径；`run_meta.json` 保存文件哈希、版本与实际样本数。程序不使用报告里写死的“64天”。

## 三、主要结果表如何读

### 1. 光谱表

`01_spectral_contrasts.csv` 中：

- `original_all`：原分区主模型，本次105天。
- `original_strict`：原分区严格模型，本次48天。
- `refined_all` / `refined_strict`：拆分过渡带，保留其余区域的补充模型。
- `seaward_minus_core`：向海林缘相对核心区的额外潮位斜率。
- `transition_minus_core`：原过渡带相对核心区的额外潮位斜率。
- `seaward_minus_transition`：两种响应的**直接差异检验**；不是比较两个区域各自是否显著。
- `transition_seaward_minus_core`：靠潮汐区一侧过渡带相对核心区的额外斜率。
- `seaward_minus_transition_seaward`：向海林缘与上述过渡带的直接差异。

单位是指数值或反射率值/米，不是百分比。模型仍有点固定效应、日期固定效应、分区×月份控制，不能估计全林共同的绝对光谱潮位斜率。

补充推断采用日期×空间区块双向聚类协方差，t检验及t区间的自由度为较小簇数减1；**未重新运行原wild bootstrap**。因此系数及SE可复现旧主模型，但本表置信区间口径和旧报告略有差别，新p值也不能标成WCR。

`q_fdr_six_indices` 在每个模型、每个对比的6个指标内校正。若某行 `variance_status=negative_two_way_variance`，其有限样本双向聚类方差为负，该行不报告SE、区间或显著性，不把负方差截为零。目前原严格模型的NDVI“向海−过渡”对比存在这一情况；不要据其系数自行判断显著。

### 2. 高低潮表

`02_composite_local_whole.csv` 中：

| 列 | 含义 |
|---|---|
| `low`、`high` | 该区域低、高潮情景下的输出 |
| `delta` | 高潮减低潮，单位见 `unit` |
| `relative_change_pct` | `100×delta/该区域low`；比较该区域自身的敏感程度 |
| `contribution_to_whole_change_pp` | `100×delta/全林low`；衡量该区域对全林总变化的贡献 |

例如，各区域的 `relative_change_pct` 不能相加；互不重叠区域的 `contribution_to_whole_change_pp` 可以相加，得到全林的变化百分比。

`part=focal` 为向海林缘＋过渡带；`rest` 为其余森林。`focal+rest=whole`，不能再加上focal内部的分区，以免重复计算。图表中中文名称见 `part_cn`。

`metric`：

- `mapped_area_ha`：识别面积，单位ha。
- `biomass_fixed_mask`：固定v4空间范围内的地上生物量汇总，不等于固定模型。RF固定/重训练分别见method。
- `biomass_in_mapped_area`：在各制图结果范围内汇总地上生物量。
- `aboveground_carbon_in_mapped_area`：上一项实际乘0.47，单位Mg C。它不是生态系统总碳，不能包含土壤碳。

仅当参考字段为 `agb_Mg_per_ha` 才输出碳量；若原GEE退回冠层高度，代码不会把高度当生物量换算为碳。

`rf_density_retrained` 指高潮影像重训练RF回归，空间统计范围仍固定。`RF_retrained_classifier` 碳量项仅重训练分类器，密度仍由低潮训练的RF回归预测，沿用原g08设计。

`03_contribution_checks.csv` 中的 `cancellation_fraction` 计算为：`1−|各分区净变化之和|/各分区变化绝对值之和`。0表示同向或没有抵消；越接近1表示正负变化抵消越多。它是输出变化的算术描述，不是潮位因果贡献率。

### 3. 逐日面积与前后期配对

逐日面积沿用原g09/s16公式：每类有效像元中判为红树林的比例，乘该类固定参考面积；整体有效覆盖至少80%，且核心区R1比例至少80%。这不保证每个局部分区都覆盖充分，所以代码另给出向海林缘和过渡带各自覆盖至少80%的敏感性分析。

`05_daily_area_effects.csv`：

- `model=original_qc`：复用旧质控，本次72天。
- `model=focal_coverage_ge80pct`：再要求局部覆盖率，本次43天；这43天**不是实测潮位版的43天严格光谱日期**。
- `slope_part_pp_per_m`：每米潮位对应的面积变化，占该分区固定参考面积的百分点。
- `slope_ha_per_m`：同一斜率换算成ha/m。
- `contribution_whole_pp_per_m`：换成统一全林参考面积分母后，该分区的斜率贡献。

回归控制年份和月份，沿用HC1标准误；这不能完全解决时间自相关或未观测混杂。年份效应不等于已验证的真实年际森林变化。

前后期配对：前期≤2020年，后期≥2023年，月份相差≤1（含12月/1月）。两期差均为后期减前期。

`06_pair_summary.csv` 中：

- `mean_part_reference_pp`：面积差除以该分区固定参考面积。
- `mean_whole_reference_pp`：面积差除以全林固定参考面积。
- `n_pairs` 是配对数，不是独立日期数；同一天会重复参加多个配对，故本表不做普通独立样本t检验。
- 分箱为右闭区间：`(-∞,-1]`、`(-1,-0.3]`、`(-0.3,0.3]`、`(0.3,1]`、`(1,∞)`；以该边界定义为准。

这些是表观变化，不能全称为“虚假损失”。

## 四、若要真正聚焦向海一侧过渡带，额外跑这一轮GEE

现有汇总CSV把所有过渡带加在一起，无法在Origin或SPSS里凭汇总值拆开。样点光谱可以直接细分，**整个区域的细分面积与碳量需要重新按新分区导出**。

### 1. 生成两个新脚本

**针对这次上传的FES结果，也提供了可直接复制到GEE的两个脚本：`gee_fes105/`。** 它们使用105个主日期、30个低潮训练日期及原项目资产路径。已核对本次上传的光谱、g08汇总和g09逐日原始CSV与对应原项目输入的SHA256完全一致。这是本次FES版本的配套脚本，不用于实测99天版本或其他站点。你若仍使用原资产路径，可直接运行这两个JS；若资产已迁移，先核对ASSET_ROOT。

如果希望从你自己的对应版本脚本生成，再使用下面的命令：

确保 `sites/cangio/gee/` 中有已填入真实资产路径和日期的原g08、g09脚本。运行：

```bat
python supplement_local_whole/make_supplement_gee.py --site-dir sites/cangio
```

会生成：

```text
sites/cangio/gee/cg_g08_tide_composite_experiment_local_whole.js
sites/cangio/gee/cg_g09_change_detection_sim_local_whole.js
```

原脚本保持原样；新脚本只修改统计分组、输出名称并增加共同范围导出。原日期、特征、阈值、训练标签、训练点和模型参数保持原设计。编码1改指靠潮汐区过渡带，编码5为其余过渡带，其他编码含义不变。不必重新运行g01—g07或生成新的分区资产。

生成器还会核对g08/g09低潮训练日期是否一致、g09主日期是否与已有原始逐日导出一致。若本地gee目录后来被实测99天版本更新，而inputs仍是FES105天，则会报版本不匹配；不要绕过检查或仅改文件名，改用上述配套FES脚本或恢复同版本源脚本。

如果 `gee/` 下同时存了多个版本而程序提示不唯一，可明确指定：

```bat
python supplement_local_whole/make_supplement_gee.py --site-dir sites/cangio --g08 sites/cangio/gee/cg_g08_tide_composite_experiment.js --g09 sites/cangio/gee/cg_g09_change_detection_sim.js
```

### 2. GEE中运行并下载新文件

分别复制新脚本到GEE Code Editor运行，在Tasks中提交导出。下载到 `sites/cangio/inputs/`，保留以下新名字：

```text
cg_s15_summary_local_whole.csv
cg_s16_daily_maps_local_whole.csv
cg_s17_factorial_common.csv
```

新像元样本 `cg_s15_pixels_local_whole.csv` 可保留，当前Python补充不依赖它。不要覆盖旧导出，也不要把旧CSV重命名伪装成新分区结果。

### 3. 汇总细分分区结果

```bat
python supplement_local_whole/local_whole.py --site-dir sites/cangio --refined-exports --skip-spectra
```

输出到 `results/6_local_whole_refined/`。光谱拆分结果已经在第一轮生成，所以这里可以跳过重复拟合。此时focal真正表示“向海林缘＋靠潮汐区一侧过渡带”，rest包含其余过渡带。全林仍包含全部森林。

### 4. 新增的碳量分解回答什么

新g08在低高潮影像**共同有效范围**内，固定低潮训练的RF回归，交叉组合：低潮范围×低潮密度、高潮范围×低潮密度、低潮范围×高潮密度、高潮范围×高潮密度。

输出 `07_carbon_factorial_decomposition.csv`：

- `mask_component`：识别范围变化对应的对称分解项。
- `density_component`：预测密度变化对应的对称分解项。
- 二者相加严格等于该共同范围内的总碳量差。

交互项对称分配给两项，因此这不是唯一的分解约定，也不是生态因果归因。`07_common_footprint_coverage.csv` 报告每类共同范围保留比例；先检查覆盖率，不能把共同范围结果直接说成完整全林结果。这个新口径与旧g08各自有效范围的汇总值可能不同。

## 五、Origin如何作图

建议在Origin中导入对应CSV，选择UTF-8编码。代码已经算好系数、区间和分区汇总；不要在Origin里再把样点日期记录当作独立样本重新做普通回归。

### 图A：局部光谱响应

1. 导入 `01_spectral_contrasts.csv`。
2. 先筛选 `model=original_all`，以及 `seaward_minus_core`、`transition_minus_core`。
3. 用 `index` 作分类X、`estimate`作Y、`zone_cn`区分系列；将 `error_minus`、`error_plus` 指定为Y误差。
4. 作分组点图/柱图并加误差线，画Y=0参考线。不同指标量纲和范围不同，优先分面展示，不把它们混算成综合指数。
5. 新分区另筛 `refined_all`，比较向海林缘与靠潮汐区过渡带。两者是否不同，看直接对比行，不靠两条误差线是否重叠判断。

### 图B：局部敏感程度与全林汇总

1. 导入 `02_composite_local_whole.csv`。
2. 筛 `metric=mapped_area_ha`、`method=R1/R2/RF_fixed`。
3. 选向海林缘、过渡带、全林，X为 `part_cn`，Y为 `relative_change_pct`，method为分组。
4. 生物量图另筛 `metric=biomass_fixed_mask`；碳量图另筛 `aboveground_carbon_in_mapped_area`，不把不同单位放在同一Y轴。

### 图C：局部怎样贡献到全林

筛选同一metric和method，只保留 `focal`、`rest` 两行，以 `contribution_to_whole_change_pp` 画正负堆积柱图，或以 `delta` 画瀑布图。两部分相加就是全林变化。也可画五个互不重叠的原分区，但不能把focal/whole再叠加进去。

### 图D：跨期潮位不一致的应用示例

导入 `06_pair_summary.csv`，筛选一个方法和一个区域，用 `tide_bin` 作X、`mean_part_reference_pp`作Y。p10/p90是配对分布的分位数，不是均值95%置信区间。图注写明配对非独立、结果含真实森林变化及其他观测差异。

## 六、SPSS是否需要

不需要为了这项补充再把所有分析搬到SPSS。普通相关、独立样本t检验或默认OLS，不能直接替代点/日期固定效应和日期×空间区块聚类推断；也不能把174个重复使用日期的配对当作174个独立样本。

SPSS可用于查看汇总和描述分布；正式面板统计用本代码，图用Origin即可。若后续希望在SPSS另建混合模型，需要另行说明随机结构和估计目标，不能声称它与当前固定效应结果完全等价。

## 七、还缺什么实验

现有数据和可选GEE导出可以支撑“敏感区域在哪里、估计结果变化多少、局部怎样贡献到全林”的分析。若希望进一步声称提高了制图精度、消除了虚假森林损失、测到了真实生物量变化，还需：

- 独立参考样本及空间/日期独立验证，不能把训练所用v4标签再当独立真值。
- 独立确认的稳定区域，用来评估潮位不匹配产生的变化误报。
- 实测潮位版与未移时FES版在共同日期、对应分组和完整流程下比较；本站结果不能直接推广为全球适用。

## 八、验证范围

```bat
python -m unittest discover -s supplement_local_whole -p test_supplement.py -v
```

本包已在本次FES归档上实际运行，原光谱系数/标准误、原72天面积、全林面积斜率、174个配对及分箱均值复现一致。新GEE脚本通过JavaScript语法检查，但未在你的GEE账号执行；细分全区统计和共同范围分解的真实数据结果，需要你提交上述导出后才能得到。详见 `VALIDATION.md`。
