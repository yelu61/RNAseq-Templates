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
| `BATCH_VECTOR` | character vector/NULL | `NULL` | 批次向量，用于批次效应诊断 |
| `PAIR_ID` | character vector/NULL | `NULL` | 配对 ID，设置后自动改为配对设计公式 |
| `SAMPLE_EXCLUDE` | character vector | `character(0)` | 要排除的样本名 |
| `COMPARISONS` | list of 3-element vectors | — | 每个比较：`c(name, treatment, control)` |
| `THRESHOLD_GRID` | data.frame | strict/standard/loose | 多阈值 DEG 网格，必须含 `name/p_cutoff/log2fc` |
| `DEFAULT_THRESHOLD` | character | `"standard"` | 默认阈值名，必须是 THRESHOLD_GRID$name 之一 |
| `DEG_PVALUE_COLUMN` | character | `"padj"` | 用于 DEG 判断的 P 值列 |
| `DEG_LFC_COLUMN` | character | `"log2FoldChange_shrunken"` | 用于 DEG 判断的 LFC 列 |
| `GSEA_RANK_COLUMN` | character | `"stat"` | GSEA 排序列 |
| `MIN_COUNT` | numeric | `10` | 低表达过滤阈值 |
| `DESIGN_FORMULA` | formula | `~ condition` | DESeq2 设计公式；配对时改为 `~ PAIR_ID + condition` |
| `RUN_TF_ANALYSIS` | logical | `FALSE` | 是否运行 DoRothEA/VIPER TF 活性分析 |
| `RUN_COMPARECLUSTER` | logical | `TRUE` | 是否运行 compareCluster（≥3 组才有效） |
| `EXPORT_EXCEL` | logical | `TRUE` | 是否输出 `1-DEG/DEG_results.xlsx` |

## `RNAseq_limma_voom_Template.ipynb`

| 参数名 | 类型 | 默认值 | 说明 |
|--------|------|--------|------|
| `INPUT_FILE` / `INPUT_FORMAT` | — | — | 同 General |
| `SAMPLE_NAMES` / `GROUPS` / `GROUP_LEVELS` | — | — | 同 General |
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
| `LOCAL_COUNTS_FILE` / `LOCAL_TPM_FILE` / `LOCAL_CLINICAL_FILE` | character | — | 本地文件模式 |
| `GENE_ID_MAP_FILE` | character/NULL | `NULL` | TCGA ENSEMBL→symbol 映射文件，NULL 自动从 SE 提取 |
| `GENES_FOR_SURVIVAL` | character vector | `c("ICAM1")` | 做 KM/Cox 的基因 |
| `CLINICAL_VARS_FOR_KM` | character vector | `c("ajcc_pathologic_stage")` | 临床变量 KM，连续变量自动二分组 |
| `TIME_UNIT` | character | `"month"` | 生存时间单位 |

## `RNAseq_TME_Deconvolution_Template.ipynb`

| 参数名 | 类型 | 默认值 | 说明 |
|--------|------|--------|------|
| `EXPR_FILE` | character | `"./1-DEG/vsd_matrix.csv"` | VST 表达矩阵 |
| `META_FILE` | character | `"./1-DEG/colData.csv"` | 元数据 |
| `GROUP_COLUMN` | character | `"condition"` | 分组列 |
| `EXPR_IS_LOG` | logical | `TRUE` | 表达矩阵是否为 log2(TPM+1)/VST 等 log 尺度；TME 工具需要非 log 输入，会自动还原 |
| `GENE_COLUMN` | character/NULL | `NULL` | 基因名列名；NULL 使用行名/第一列 |
| `SAMPLE_COLUMN` | character | `"sample"` | 样本名列 |
| `GROUP_LEVELS` | character vector/NULL | `NULL` | 分组水平顺序 |
| `SPECIES` | character | `"human"` | `"human"` 或 `"mouse"`；小鼠数据会经 `biomaRt::getLDS` 把 MGI symbol 转为 HGNC symbol 后再进行反卷积 |
| `GROUP_COLORS` | named character vector/NULL | `NULL` | 自定义分组颜色，例如 `c("Control"="#6F6F6F", "Treatment"="#E07B54")`；NULL 时自动生成 |
| `RUN_ESTIMATE` | logical | `TRUE` | 是否运行 native ESTIMATE |
| `RUN_IOBR` | logical | `TRUE` | 是否运行 IOBR |
| `IOBR_METHODS` | character vector | `c("estimate", "cibersort", "epic", "xcell")` | 要运行的 IOBR 方法 |
| `IOBR_PERM` | numeric | `1000` | CIBERSORT permutation 次数 |
| `IOBR_ARRAYS` | logical | `FALSE` | RNA-seq 设为 FALSE；microarray 设为 TRUE |
| `RUN_CIBERSORT` | logical | `FALSE` | 是否运行 native CIBERSORT（需自备 `CIBERSORT.R` 和 `LM22.txt`） |
| `CIBERSORT_SCRIPT` | character | `"./CIBERSORT.R"` | native CIBERSORT 脚本路径 |
| `CIBERSORT_SIGNATURE` | character | `"./LM22.txt"` | CIBERSORT signature 文件路径 |

### 各 TME 方法输入要求

| 方法 | 需要的输入 | 备注 |
|------|-----------|------|
| ESTIMATE (native) | 非 log 归一化表达；行名 = HGNC symbol | 与人类基质/免疫 signature 取交集 |
| CIBERSORT (native) | 非 log 归一化表达；第一列 = HGNC symbol | `LM22.txt` 为人类 reference |
| IOBR `estimate` | 非 log TPM-like；行名 = HGNC symbol | 包装自 ESTIMATE |
| IOBR `cibersort` | 非 log TPM-like；行名 = HGNC symbol | RNA-seq 设 `arrays = FALSE` |
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
