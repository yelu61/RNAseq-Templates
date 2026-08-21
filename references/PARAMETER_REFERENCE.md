# 参数速查表 (Parameter Reference)

本文件汇总了各 notebook 中常用的配置参数。每个参数包含名称、类型、默认值和用途说明。

## `RNAseq_General.ipynb`

| 参数名 | 类型 | 默认值 | 说明 |
|--------|------|--------|------|
| `SPECIES` | character | `"mouse"` | `"mouse"` 或 `"human"`，决定 org.Mm.eg.db / org.Hs.eg.db |
| `INPUT_FILE` | character | `"./0-Data/...tsv"` | 原始 count 矩阵路径 |
| `INPUT_FORMAT` | character | `"tsv"` | `"tsv"` / `"csv"` / `"excel"` |
| `GENE_NAME_COL` | character | `"gene_name"` | 基因名列名 |
| `BIOTYPE_COL` | character/NULL | `"gene_biotype"` | 基因类型列名，NULL 表示不过滤 |
| `BIOTYPE_FILTER` | character | `"protein_coding"` | 保留的基因类型 |
| `SAMPLE_NAMES` | character vector | — | 样本名，必须与 count 列名一致 |
| `GROUPS` | character vector | — | 每个样本所属分组 |
| `GROUP_LEVELS` | character vector | — | 分组水平（控制顺序和参考组） |
| `BATCH_VECTOR` | character vector/NULL | `NULL` | 批次向量，用于批次效应诊断；若提供但 `DESIGN_FORMULA` 未含 `batch` 会 warning |
| `PAIR_ID` | character vector/NULL | `NULL` | 配对 ID，设置后自动改为配对设计公式 |
| `SAMPLE_EXCLUDE` | character vector | `character(0)` | 要排除的样本名 |
| `COMPARISONS` | list of 3-element vectors | — | 每个比较：`c(name, treatment, control)` |
| `THRESHOLD_GRID` | data.frame | strict/standard/loose | 多阈值 DEG 网格，必须含 `name/p_cutoff/log2fc`；可选 `p_column` 按层指定 `padj` 或仅探索用的 `pvalue` |
| `DEFAULT_THRESHOLD` | character | `"standard"` | 默认阈值名，必须是 THRESHOLD_GRID$name 之一 |
| `DEG_PVALUE_COLUMN` | character | `"padj"` | 用于 DEG 判断的 P 值列 |
| `DEFAULT_DEG_PVALUE_COLUMN` | character | 由默认阈值解析 | 默认层实际使用的 P 值列；由该行 `p_column` 覆盖，否则继承 `DEG_PVALUE_COLUMN` |
| `DEG_LFC_COLUMN` | character | `"log2FoldChange_raw"` | 用于 DEG 判断的 LFC 列；可改为 shrunken 做更保守展示 |
| `GSEA_RANK_COLUMN` | character | `"stat"` | GSEA 排序列 |
| `PAIRWISE_P_ADJUST_METHOD` | character | `"BH"` | GSVA/单基因两两比较图的 P 值校正方法 |
| `MIN_COUNT` | numeric | `10` | 低表达过滤阈值 |
| `DESIGN_FORMULA` | formula | `~ condition` | DESeq2 设计公式；配对时改为 `~ PAIR_ID + condition` |
| `RUN_TF_ANALYSIS` | logical | `FALSE` | 是否运行 DoRothEA/VIPER TF 活性分析 |
| `RUN_COMPARECLUSTER` | logical | `TRUE` | 是否运行 compareCluster（≥3 组才有效） |
| `EXPORT_EXCEL` | logical | `TRUE` | 是否输出 `1-DEG/DEG_results.xlsx` |
| `GENERATE_HTML_REPORT` | logical | `TRUE` | 是否渲染统一 HTML 报告 `RNAseq_report.html`（无 quarto CLI 时自动回退 rmarkdown 渲染 `.Rmd` 双胞胎模板） |
| `REPORT_TITLE` | character | `"RNA-seq Analysis Report"` | HTML 报告标题 |
| `SCAFFOLD_REPORT_INTERPRETATION` | logical | `FALSE` | 首次需要项目说明时设为 TRUE，创建不覆盖的 `report_interpretation/<section>.md` 注释骨架 |

`p_column = "pvalue"` 的层只用于假设生成，不能作为项目主结论阈值；相应 ORA 必须由 adjusted-P 层或 GSEA/GSVA 等等级证据确认。

## `RNAseq_limma_voom_Template.ipynb`

| 参数名 | 类型 | 默认值 | 说明 |
|--------|------|--------|------|
| `INPUT_FILE` / `INPUT_FORMAT` | — | — | 同 General |
| `SAMPLE_NAMES` / `GROUPS` / `GROUP_LEVELS` | — | — | 同 General |
| `SPECIES` | character | `"human"` | `"human"` 或 `"mouse"`；决定富集分析使用的 org.db 和 KEGG organism code |
| `BATCH_VECTOR` | character vector/NULL | `NULL` | 批次向量，传给 `limma::removeBatchEffect` |
| `COMPARISONS` | list | — | 同 General |
| `DEG_PADJ_CUTOFF` | numeric | `0.05` | DEG 阈值 |
| `DEG_LFC_CUTOFF` | numeric | `0.5` | DEG log2FC 阈值 |
| `MIN_COUNT` | numeric | `10` | filterByExpr 最小 count |
| `MIN_SAMPLE_FRAC` | numeric | `0.5` | filterByExpr 最小样本比例 |

## `RNAseq_TimeCourse_Template.ipynb`

| 参数名 | 类型 | 默认值 | 说明 |
|--------|------|--------|------|
| `EXPR_FILE` | character | `"./1-DEG/vsd_matrix.csv"` | VST 表达矩阵，用于 Mfuzz |
| `META_FILE` | character | `"./1-DEG/colData.csv"` | 样本元数据，含 time / condition |
| `GENE_COLUMN` | character/NULL | `NULL` | 基因名列名；NULL 使用行名/第一列 |
| `SAMPLE_COLUMN` | character | `"sample"` | 样本名列 |
| `TIME_COLUMN` | character | `"time"` | 时间列名 |
| `GROUP_COLUMN` | character/NULL | `"condition"` | 分组列名，可选 |
| `TIME_LEVELS` | character vector/NULL | `NULL` | 时间水平顺序 |
| `RAW_COUNTS_FILE` | character | `"./0-Data/raw_counts.tsv"` | 原始 count 矩阵，用于 Section 6 |
| `COUNT_META_FILE` | character | `"./0-Data/metadata.csv"` | 原始 count 的样本元数据 |
| `SUBJECT_COL` | character/NULL | `NULL` | 受试者 ID，用于配对设计 |
| `RUN_TIMEPOINT_DEG` | logical | `TRUE` | 是否运行时间点 vs baseline DEG |
| `BASELINE_TIME` | character/NULL | `NULL` | 基线时间点，NULL 取 TIME_LEVELS 第一个 |
| `MFUZZ_N_CLUSTERS` | numeric | `5` | Mfuzz 聚类数 |
| `MFUZZ_MIN_ACORE` | numeric | `0.7` | Mfuzz core gene 阈值 |

## `RNAseq_TCGA_GEO_Template.ipynb`

| 参数名 | 类型 | 默认值 | 说明 |
|--------|------|--------|------|
| `DOWNLOAD_FROM_GDC` | logical | `TRUE` | 是否从 TCGA GDC 下载 |
| `DOWNLOAD_FROM_GEO` | logical | `FALSE` | 是否从 GEO 下载 SeriesMatrix |
| `GEO_ACCESSION` | character | `"GSE12345"` | GEO accession |
| `TCGA_PROJECT` | character | `"TCGA-STAD"` | TCGA project ID |
| `GDC_COUNTS_ASSAY` | character/NULL | `NULL` | GDC SummarizedExperiment 中 count assay 名称；NULL 自动检测 |
| `GDC_TPM_ASSAY` | character/NULL | `NULL` | GDC SummarizedExperiment 中 TPM assay 名称；NULL 自动检测 |
| `LOCAL_COUNTS_FILE` / `LOCAL_TPM_FILE` / `LOCAL_CLINICAL_FILE` | character | — | 本地文件模式 |
| `LOCAL_GENE_COLUMN` | character/NULL | `NULL` | 本地矩阵的基因 ID 列；数值型 Entrez ID 必须显式填写 |
| `GENE_ID_MAP_FILE` | character/NULL | `NULL` | TCGA ENSEMBL→symbol 映射文件，NULL 自动从 SE 提取 |
| `GENES_FOR_SURVIVAL` | character vector | `c("ICAM1")` | 做 KM/Cox 的基因 |
| `CLINICAL_VARS_FOR_KM` | character vector | `c("ajcc_pathologic_stage")` | 临床变量 KM，连续变量自动二分组 |
| `TIME_UNIT` | character | `"month"` | 生存时间单位 |

## `RNAseq_TME_Deconvolution_Template.ipynb`

| 参数名 | 类型 | 默认值 | 说明 |
|--------|------|--------|------|
| `INPUT_MODE` | character | `"raw_counts"` | 推荐默认；`expression` 仅用于已有 TPM/log2(TPM+1) |
| `RAW_COUNTS_FILE` | character | `"./0-Data/featureCounts_merged_count.annot.tsv"` | 原始整数 count 表；默认输入 |
| `RAW_COUNTS_FORMAT` | character | `"tsv"` | `tsv`、`csv` 或 `excel` |
| `EXPR_FILE` | character | `"./0-Data/TPM_matrix.csv"` | 仅在 expression 模式使用；不能提供 VST/rlog |
| `EXPR_UNIT` | character | `"tpm"` | `tpm` 或 `log2_tpm`；必须显式声明 |
| `META_FILE` | character | `"./1-DEG/colData.csv"` | 元数据 |
| `GENE_COLUMN` | character/NULL | `"gene_id"` | TPM 计算优先使用唯一稳定 ID；Ensembl ID 后续转 symbol；数值型 ID 必须显式指定列名 |
| `GENE_LENGTH_COLUMN` | character/NULL | `NULL` | raw count 表中基因长度列名（bp 或 kb，由 `GENE_LENGTH_UNIT` 控制） |
| `GENE_LENGTH_UNIT` | character | `"bp"` | `"bp"` 或 `"kb"` |
| `GENE_START_COL` / `GENE_END_COL` | character | `"gene_start"` / `"gene_end"` | 若未提供 `GENE_LENGTH_COLUMN`，则用起止坐标计算长度 |
| `SAMPLE_COLUMN` | character | `"sample"` | 样本名列 |
| `GROUP_COLUMN` | character | `"condition"` | 分组列 |
| `GROUP_LEVELS` | character vector/NULL | `NULL` | 分组水平顺序 |
| `SPECIES` | character | `"human"` | `"human"` 或 `"mouse"`；小鼠数据会经 `convert_expression_rownames()` 把 MGI symbol 转为 HGNC symbol 后再进行反卷积 |
| `GROUP_COLORS` | named character vector/NULL | `NULL` | 自定义分组颜色；NULL 时自动生成 |
| `RUN_ESTIMATE` | logical | `TRUE` | 是否运行 native ESTIMATE |
| `RUN_IOBR` | logical | `FALSE` | 是否运行 IOBR；首次运行可能下载方法参考数据，离线生产运行前需先完成缓存预检 |
| `IOBR_METHODS` | character vector | `c("estimate", "cibersort", "epic", "xcell")` | 要运行的 IOBR 方法 |
| `IOBR_PERM` | numeric | `1000` | CIBERSORT permutation 次数 |
| `IOBR_ARRAYS` | logical | `FALSE` | RNA-seq 设为 FALSE；microarray 设为 TRUE |
| `RUN_CIBERSORT` | logical | `FALSE` | 是否运行 native CIBERSORT（项目已内置 `references/CIBERSORT/` 资源） |
| `CIBERSORT_SCRIPT` | character | `NULL`（自动定位） | native CIBERSORT 脚本路径；`NULL` 时从项目 git 根目录自动定位 `references/CIBERSORT/CIBERSORT.R` |
| `CIBERSORT_SIGNATURE` | character | `NULL`（自动定位） | CIBERSORT signature 文件路径；`NULL` 时按 `SPECIES` 自动选择 `LM22.txt`（human）或 `cibersort_mouse_22.csv`（mouse） |
| `CIBERSORT_PERM` | numeric | `1000` | CIBERSORT permutation 次数（native 与 IOBR 一致） |
| `CIBERSORT_QN` | logical | `FALSE` | RNA-seq 关闭分位数归一化；microarray 通常设为 TRUE |
| `RUN_CIBERSORT_COMPARISON` | logical | `TRUE` | 当 `RUN_CIBERSORT` 与 `RUN_IOBR` 均含 `"cibersort"` 时，是否自动生成 native vs IOBR 对比表和 PDF（输出到 `OUTDIR`） |

### 各 TME 方法输入要求

| 方法 | 需要的输入 | 备注 |
|------|-----------|------|
| ESTIMATE (native) | 非 log 归一化表达；行名 = HGNC symbol | 与人类基质/免疫 signature 取交集 |
| CIBERSORT (native) | 非 log 归一化表达；行名 = HGNC symbol（human）或 MGI symbol（mouse） | 内置 `LM22.txt` 为人类 reference；`cibersort_mouse_22.csv` 为小鼠 reference，可直接配合小鼠基因符号使用 |
| IOBR `estimate` | 非 log TPM-like；行名 = HGNC symbol | 包装自 ESTIMATE |
| IOBR `cibersort` | 非 log TPM-like；行名 = HGNC symbol | RNA-seq 设 `arrays = FALSE`；小鼠数据需先经 `prepare_tme_expression(species="mouse")` 转换为 HGNC symbol |
| IOBR `epic` | 非 log TPM-like；行名 = HGNC symbol | 输出细胞比例，总和 ≤ 1 |
| IOBR `xcell` | 非 log TPM-like；行名 = HGNC symbol | 输出富集分数，非比例 |
| ssGSEA | log 或非 log 均可；行名与 signature 匹配 | 内置免疫 signature 为人类基因符号 |

## `RNAseq_WGCNA_Template.ipynb`

| 参数名 | 类型 | 默认值 | 说明 |
|--------|------|--------|------|
| `EXPR_FILE` / `META_FILE` | — | — | 同 TME |
| `SOFT_THRESHOLD_POWERS` | numeric vector | `c(1:20)` | 软阈值候选幂 |
| `MIN_MODULE_SIZE` | numeric | `30` | 最小模块大小 |
| `MERGE_CUT_HEIGHT` | numeric | `0.25` | 模块合并高度 |
| `TRAIT_COLUMNS` | character vector | — | 与模块相关的性状/临床变量列 |

## 通用注意事项

1. **样本名必须完全匹配**：`SAMPLE_NAMES` 或 metadata 中的 sample 列必须与 count 矩阵列名一致（包括大小写）。
2. **GROUP_LEVELS 顺序**：第一个水平通常作为对照组/参考组，影响 DESeq2 / limma 的结果方向。
3. **比较方向**：每个 `COMPARISONS` 条目格式为 `c("name", "treatment", "control")`，结果表示 treatment vs control。
4. **Excel 导出**：依赖 `openxlsx`，已通过 `install_dependencies.R` 安装。

## General runner 生命周期参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `RUN_ROLE` | `"candidate"` | `candidate`、`canonical`、`sensitivity`、`repro_check` 或 `superseded`；不改变统计模型 |
| `PARENT_RUN_ID` | `NA_character_` | sensitivity/replacement run 的来源 run_id |
| `RUN_CHANGE_NOTE` | `""` | 与 parent 相比的唯一、简短变更原因 |
| `RUN_RETENTION` | `"full"` | `full`、`slim` 或 `metadata_only`；模板仅记录，不自动清理 |

General runner 在完成时写 `run_manifest.csv`。使用
`Rscript tools/build_run_registry.R analysis/runs` 汇总并识别完全重复 run。
