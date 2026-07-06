# 零基础上手指南 (Getting Started)

如果你会操作 Excel，就会使用这个 RNA-seq 分析模板。本指南面向**完全没有编程经验**的用户。

---

## 你需要准备的

1. 一台电脑（Windows / macOS / Linux 都可以）
2. 已经安装好的 **R** 和 **RStudio**（见下方安装步骤）
3. 你的 RNA-seq count 表达矩阵文件（基因 × 样本）
4. 你的样本分组信息

---

## 第一步：安装 R 和 RStudio

### 1.1 安装 R

1. 打开浏览器，访问：https://cloud.r-project.org/
2. 点击与你电脑系统对应的链接：
   - Windows 用户：Download R for Windows → base → Download R-4.x.x for Windows
   - macOS 用户：Download R for macOS → 选择第一个 `.pkg` 文件下载安装
3. 下载后双击安装，全部点"继续"即可

### 1.2 安装 RStudio

1. 访问：https://posit.co/download/rstudio-desktop/
2. 点击 **Download RStudio Desktop**
3. 下载后双击安装

安装完成后，打开 RStudio，你应该看到一个有 4 个面板的窗口。

---

## 第二步：安装分析所需的所有 R 包

1. 下载或克隆本仓库到电脑本地
2. 在 RStudio 菜单点击 **File → Open Project...**，选择仓库文件夹
3. 在 RStudio 左下角找到 **Console** 面板（有一个 `>` 提示符）
4. 复制粘贴下面这行代码，按回车：

```r
Rscript install_dependencies.R
```

如果提示找不到 `Rscript`：
- 在 RStudio Console 中输入：

```r
source("install_dependencies.R")
```

这会开始下载安装所有需要的包，**第一次可能需要 10-30 分钟**，请耐心等待，不要关闭 RStudio。

---

## 第三步：运行演示（Demo）

仓库里已经准备好了可以一键跑通的演示数据。

1. 在 RStudio 右下角 **Files** 面板中，依次打开：
   ```
   examples → demo_RNAseq_General → RNAseq_General.ipynb
   ```
2. 如果提示需要安装 Jupyter/IRkernel，按照提示安装即可
3. 打开后，点击菜单 **Cell → Run All**
4. 等待几分钟，左侧会生成三个文件夹：
   - `1-DEG/`：差异表达结果
   - `2-GSEA/`：富集分析结果
   - `3-Visualization/`：所有 PDF 图片

如果演示能跑通，说明环境没有问题。

---

## 第四步：换成你自己的数据

### 4.1 准备输入文件

#### 文件 1：表达矩阵（必需）

格式示例（`counts.tsv`）：

| gene_name | gene_biotype | Sample_A | Sample_B | Sample_C | Sample_D |
|-----------|--------------|----------|----------|----------|----------|
| Gene1     | protein_coding | 100      | 120      | 80       | 200      |
| Gene2     | protein_coding | 50       | 45       | 60       | 55       |

要求：
- 必须是**原始 count 值**（整数，不能是 TPM/FPKM）
- 必须有一列基因名（如 `gene_name`）
- 如果有 `gene_biotype` 列，可以自动过滤蛋白编码基因
- 文件可以是 `.tsv`、`.csv` 或 `.xlsx`

#### 文件 2：样本分组信息（可选但推荐）

格式示例（`metadata.csv`）：

| sample    | condition |
|-----------|-----------|
| Sample_A  | Control   |
| Sample_B  | Control   |
| Sample_C  | Treatment |
| Sample_D  | Treatment |

### 4.2 复制演示文件夹

1. 在 RStudio 右下角 **Files** 面板中，右键点击 `examples/demo_RNAseq_General/`
2. 选择 **Copy**，然后粘贴到你的项目文件夹（比如桌面新建一个文件夹）
3. 把刚才准备好的表达矩阵和 metadata 文件放进 `0-Data/` 文件夹

### 4.3 修改参数

1. 打开复制出来的 `RNAseq_General.ipynb`
2. 找到第 2 个单元格 **1. Parameter Configuration**
3. 只需要修改以下几行：

```r
SPECIES       <- "mouse"                    # 改成 "mouse" 或 "human"
INPUT_FILE    <- "./0-Data/你的counts.tsv"   # 改成你的文件名
INPUT_FORMAT  <- "tsv"                      # tsv / csv / excel
GENE_NAME_COL <- "gene_name"                # 你的基因名列名
BIOTYPE_COL   <- "gene_biotype"             # 如果没有这一列，改成 NULL

SAMPLE_NAMES <- c("Sample_A", "Sample_B", "Sample_C", "Sample_D")
GROUPS       <- c("Control", "Control", "Treatment", "Treatment")
GROUP_LEVELS <- c("Control", "Treatment")

COMPARISONS <- list(
  c("Treatment_vs_Control", "Treatment", "Control")
)
```

### 4.4 运行分析

点击菜单 **Cell → Run All**，等待完成即可。

---

## 第五步：查看结果

跑完后，你会在当前文件夹看到：

```text
1-DEG/                          <- 差异表达表格
2-GSEA/                         <- 富集分析结果
3-Visualization/                <- PDF 图片
Analysis_summary.txt            <- 文字版分析摘要
sessionInfo.txt                 <- 运行环境记录
```

### 重要输出文件

| 文件/文件夹 | 含义 |
|-------------|------|
| `1-DEG/DEG_threshold_summary.csv` | 每个比较上下调基因数量 |
| `1-DEG/all_genes/` | 所有基因的 DESeq2 完整结果 |
| `3-Visualization/Volcano_*.pdf` | 火山图 |
| `3-Visualization/DEG_heatmap.pdf` | 差异基因热图 |
| `3-Visualization/*_ORA_*.pdf` | GO/KEGG 富集图 |
| `3-Visualization/*_GSEA_*.pdf` | GSEA 富集图 |

---

## 常见问题

### Q1: 运行中报错"找不到文件"

检查 `INPUT_FILE` 的路径是否正确。如果是放在 `0-Data/` 文件夹里，应该写成：

```r
INPUT_FILE <- "./0-Data/your_file.tsv"
```

### Q2: 我不知道分组怎么写

假设你有 3 个对照和 3 个处理：

```r
SAMPLE_NAMES <- c("Ctrl1", "Ctrl2", "Ctrl3", "Treat1", "Treat2", "Treat3")
GROUPS       <- c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment")
GROUP_LEVELS <- c("Control", "Treatment")
```

`SAMPLE_NAMES` 必须与表达矩阵里的列名完全一致。

### Q3: 我只有一个比较，不需要 compareCluster

把下面这行改成 `FALSE`：

```r
RUN_COMPARECLUSTER <- FALSE
```

### Q4: 安装依赖时某个包失败

在 RStudio Console 中单独运行：

```r
install.packages("包名")
```

如果提示 Bioconductor 包失败：

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("包名")
```

---

## 验证安装

安装完成后，建议运行仓库自带的 smoke test 确认环境和核心流程正常：

```r
Rscript examples/run_demo_smoke_test.R
```

如果显示 `Smoke test PASSED`，说明所有依赖和核心分析流程都可以正常工作。

> 注意：demo 数据使用真实小鼠基因 symbol。如果 enrichment 步骤产生的显著条目较少，某些阈值下可能没有 GO/KEGG ORA 结果文件；smoke test 会检查 ORA 汇总表是否存在，而不是要求某个特定阈值目录下一定产生结果。

---

## 下一步

跑通通用流程后，根据你的分析目的选择专题模板。更详细的对比和决策流程请参考 [references/TEMPLATE_SELECTION.md](references/TEMPLATE_SELECTION.md)。

| 分析目的 | 使用模板 |
|----------|----------|
| 两组/多组差异表达 | `RNAseq_General.ipynb` |
| 时间序列/药物处理梯度 | `RNAseq_TimeCourse_Template.ipynb` |
| 肿瘤微环境免疫细胞浸润 | `RNAseq_TME_Deconvolution_Template.ipynb` |
| 共表达网络/WGCNA | `RNAseq_WGCNA_Template.ipynb` |
| TCGA/GEO 公共数据挖掘 | `RNAseq_TCGA_GEO_Template.ipynb` |
| limma-voom 替代 DESeq2 | `RNAseq_limma_voom_Template.ipynb` |

---

## 常见问题

更完整的故障排查请查看 [references/TROUBLESHOOTING.md](references/TROUBLESHOOTING.md)。

---

## 需要帮助？

如果按本指南操作仍然无法运行，请准备以下信息寻求帮助：

1. 你操作到哪一步出错
2. 完整的错误提示（红色文字）
3. 你的参数配置截图或复制文本
4. `sessionInfo.txt` 文件内容
