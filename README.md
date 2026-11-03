# NTSminer

**一句话**：面向复杂环境样品的 **非靶向筛选（NTS）** 与 **双网络传播**（结构 + 质谱），整合 LC-MS(/MS) 处理、网络构建与可视化，用于 PFAS 及其转化产物的发现、溯源与注释的 R 包。

- **当前版本**：以根目录 [`DESCRIPTION`](DESCRIPTION) 为唯一来源；版本变更记录见 [`NEWS.md`](NEWS.md)（不在本文档重复历史版本列表）。
- **License**：MIT（见根目录 `LICENSE`）。

## 安装

在 R 会话中（需已安装 **[Bioconductor](https://bioconductor.org/)** 相关依赖及本包 `Imports` 所列包）：

```r
# 从本地源码目录安装（开发常用）
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")
remotes::install_local("路径/到/NTSminer", dependencies = TRUE)
```

或直接于克隆目录内调试加载：

```r
pkgload::load_all()   # 或 devtools::load_all()
```

## 使用自己的数据

本发布版不包含演示原始数据。请将 MS1 / MS2 路径替换为自己的数据目录；包内仍提供常用 PFAS 参考库、诊断碎片表及模型资源。

```r
ms1_dir <- "path/to/MS1"
ms2_dir <- "path/to/MS2"

# 三套常用 PFAS 谱库（in-house L1 / OECD online L2a / OECD in-silico L2b）
db <- load_builtin_pfas_library()  # 默认合并三库

# 拿任一内置文件的磁盘路径
bundled_nts_file("extdata/tables/diagnostic_fragments.csv")
```

内置表、谱库和模型通过 `bundled_nts_file()` 或相应的包内加载函数定位；MS1/MS2 原始文件由用户提供。

## 运行与检查

```r
# 启动附带 Shiny 界面
NTSminer::launch_ntsminer()

# 包检查（需开发工具链）
devtools::check()
devtools::test()
```

Vignette 入口：[vignettes/ntsminer-workflow.Rmd](vignettes/ntsminer-workflow.Rmd)。

## 环境归趋 / ADMET / ToxPi

| 步骤 | 函数 | 说明 |
|------|------|------|
| 1 | **`canonicalize_smiles()`** | 预测前统一 SMILES，避免 EPI 改写导致合并错位 |
| 2 | **`predict_episuite()`** / **`predict_admet_risk()`** | EPI Suite REST + ADMET-AI（或 ADMETlab 3.0） |
| 3 | **`prepare_toxpi()`** | 六切片得分（0–1）；传 **`compound_ids`**（如筛查表 `name` / `Lab_ID`），勿只靠 `episuite_Name` |
| 4 | **`run_toxpi()`** / **`plot_toxpi()`** | 需 **`toxpiR`**（`Suggests`） |

## logIE 预测与半定量 RIE

| 步骤 | 函数 | 说明 |
|------|------|------|
| 1 | **`predict_logIE()`** | 纯 rcdk 指纹+描述符 XGBoost 预测 logIE（POS R²=0.565 / NEG R²=0.644） |
| 2 | **`predict_rie()`** | SMILES → logIE 模型 → `rie = 10^pred_logIE`，自动加载模型/检测 SMILES 列 |
| 3 | **`run_semi_quantification()`** | 内标 → 归一化 → RIE → 浓度，透传 `logIE_model`/`polarity`/`smiles_col` |
| - | **`build_logIE_model()`** | 训练自定义 logIE 模型（XGBoost/RF），保存/加载 S3 `nts_logIE_model` |
| - | **`predict_organic_percentage()`** | RT → 有机相比例插值（`R/04_feature_filter.R`），支持 linear/constant/spline |

```r
model <- load_logIE_model(bundled_nts_file("extdata/models/logIE/logIE_model_inhouse_pos.rds"))
pred  <- predict_logIE(smiles_vec, model = model, polarity = "positive")
rie   <- predict_rie(feature_table, logIE_model = model, polarity = "positive")
```

```r
cf <- classify_classyfire(smiles_tbl$SMILES, compound_ids = smiles_tbl$name)
```

```r
devtools::load_all(".")
install.packages("toxpiR")  # 作图与打分

smiles_tbl <- utils::read.csv("path/to/compound_smiles.csv", stringsAsFactors = FALSE)
smiles_canon <- canonicalize_smiles(smiles_tbl$SMILES)
epi   <- predict_episuite(smiles_canon)
admet <- predict_admet_risk(smiles_canon, backend = "admet_ai")
prep  <- prepare_toxpi(epi, admet, query_smiles = smiles_canon, compound_ids = smiles_tbl$name)
res   <- run_toxpi(prep, weight_multiplier = 1)
plot_toxpi(sort(res)[seq_len(min(20L, length(res)))])
```

