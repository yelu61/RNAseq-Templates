# General 模板 — 命令行运行（run_analysis.R + config.R)

`notebooks/RNAseq_General.ipynb` 适合**交互式探索**（改 gene set 做 GSVA、临时加可视化、逐 cell 检查）。但在实际项目里，标准流程更适合**一键、可复现、可批量**地跑完——这就是 `run_analysis.R` + `config.R` 的用途。两者共享分析意图与 helper；生产 runner 是最终统计事实来源，产物不要求逐文件一致。

## 用法

真实项目优先采用 [完整运行示例](../../BEST_PRACTICES.md#1--reproducible-run)：
固定 backend 版本，只复制配置，在新建的 `analysis/runs/<run_id>/` 中调用
中央 runner，并显式设置 `RNASEQ_LIB_DIR`。以下复制方式适合独立运行目录；
不要在共享模板目录或已经完成的 run 中执行。

1. 把本目录两个文件复制到你的项目文件夹：
   ```
   my_project/
     config.R          # 只改这个
     run_analysis.R    # 不用改
     0-Data/           # 你的输入数据
   ```
2. 编辑 `config.R`：输入路径、样本名、分组、比较、阈值、GSVA gene set、KEY_GENES 等。
3. 运行：
   ```bash
   Rscript run_analysis.R                 # 用 runner 同目录的 config.R
   Rscript run_analysis.R my_config.R     # 指定另一个 config
   ```

跑完产出 `0-Config/ 1-DEG/ 2-GSEA/ 3-Visualization/ Analysis_summary.txt sessionInfo.txt run_manifest.csv`，以及（若 `GENERATE_HTML_REPORT <- TRUE`)`RNAseq_report.html`。

## 真实项目推荐目录

不要把 runner 放进 `notebooks/` 后直接运行。相对路径会以当前工作目录
为根，从而把全部 DEG、GSEA 和图形堆到 notebook 旁边。推荐每次运行使用
不可覆盖的原生产物目录，再把审阅后的成果整理进 `results/`：

```text
analysis/config/                 # 项目配置
analysis/scripts/                # runner 与项目专用脚本
analysis/runs/<run_id>/          # 本次完整原生产物
results/tables/                  # 精选、可交付表格
results/figures/                 # 精选 PDF/SVG/TIFF 主文件
results/reports/                 # HTML/Markdown 报告
results/report_assets/           # 可重建的网页预览图
```

每次新分析使用新的 `run_id`，记录模板版本、配置与输入校验和；不要覆盖旧
运行。完整约定见 [真实项目输出结构](../../references/PROJECT_OUTPUT_LAYOUT.md)。

`config.R` 的 `RUN_ROLE`、`PARENT_RUN_ID`、`RUN_CHANGE_NOTE` 和
`RUN_RETENTION` 记录 candidate/canonical/sensitivity 等生命周期信息。完成
多个 run 后执行 `Rscript /path/to/backend/tools/build_run_registry.R /path/to/project/analysis/runs`，
生成 `RUN_REGISTRY.csv`。新版 `run_inputs.csv` 校验已声明的主/辅助/参考文件，
`run_runtime.csv` 记录 R/包版本。只有完整 schema-v2 签名才参与候选重复分组；
General 的在线 KEGG 参考未留存，因此不满足完整条件。`duplicate_of` 不能
作为自动清理依据；详见 [冻结审计](../../references/FREEZE_AUDIT.md)。

## 可选分析开关(config.R)

| 开关 | 作用 | 依赖 |
|---|---|---|
| `RUN_TF_ANALYSIS` | DoRothEA/VIPER 转录因子活性 | dorothea, viper, limma；前两者需额外安装，默认安装脚本未包含 |
| `RUN_TME` | 肿瘤微环境去卷积（ESTIMATE / IOBR / ssGSEA)，输出到 `4-TME/` | 按子开关需要 IOBR、estimate、GSVA；小鼠默认使用本地 babelgene |
| `RUN_COMPARECLUSTER` | 多组 compareCluster 富集比较 | ≥3 组才有意义 |
| `EXPORT_EXCEL` | 导出多阈值 DEG 到 `1-DEG/DEG_results.xlsx` | openxlsx |
| `GENERATE_HTML_REPORT` | 渲染统一 HTML 报告 | quarto CLI 或 rmarkdown |

`THRESHOLD_GRID` 可为某一行添加 `p_column = "pvalue"`，作为低 DEG 情形下补充 ORA 的**探索层**；默认层的图、交集和报告会自动使用该层实际解析出的 P 值列。名义 P 层不能作为主结论阈值。

**TME 特别注意**：去卷积用 TPM，不能从 VST 还原。`RUN_TME=TRUE` 会重新读取全量原输入，保留最终纳入的样本，先用全部特征及匹配长度计算 TPM，再做 symbol/同源映射；DEG 和 biotype 过滤不再改变 TPM 分母。`TME_GENE_ID_COLUMN` 可指定唯一稳定 ID，默认优先 `gene_id/Geneid`，否则用 `GENE_NAME_COL`：
- 有长度列 → 设 `TME_GENE_LENGTH_COLUMN <- "gene_length"`;
- 优先使用与计数注释一致的外显子/feature 长度（如 featureCounts `Length`），不能把包含内含子的整个基因跨度直接当作同等长度。
- 起止坐标仅在其定义符合计数量化单位时作为备用；需记录长度来源。
- General 与专题 TME 现在共享映射流程：优先原表 symbol 注释，必要时做本物种 Ensembl→symbol；小鼠默认用本地 babelgene 转人源直系同源。IOBR 仍须验证方法参考缓存。
- 显式设置 `TME_ORTHOLOG_CACHE` 可继续使用旧 biomaRt/缓存路径；缺失映射可能联网，失败时缓存降级会警告。固定此选择及缓存版本，避免两次分析使用不同映射来源。
- General runner 默认 IOBR 方法为 `estimate/cibersort/epic/mcpcounter`；`xcell` 需要在线 biomaRt，只有在网络和参考库已验证时才显式加入。
- 子开关 `RUN_TME_ESTIMATE` / `RUN_TME_IOBR` / `RUN_TME_SSGSEA` 可分别控制三个方法。
- 导出的 `TME_gene_mapping.csv`、`TME_input_coverage.csv` 与 `ssGSEA_signature_coverage.csv` 记录转换来源和覆盖度。当前 ssGSEA 使用 `ssgseaParam()`；旧 v0.11.0 的同名结果实际为 GSVA，需新建 run 重算。General 不含本地 native CIBERSORT，需要人/鼠本地参考分支时使用专题 TME。

HTML 报告的项目说明不要写进共享模板或已生成 HTML。首次运行时设 `SCAFFOLD_REPORT_INTERPRETATION <- TRUE`，再编辑 `report_interpretation/<section>.md`；这些文件会在重渲染时保留并嵌入相应章节。用 backend 中 `tools/render_report.R` 的绝对路径独立重渲染，并显式传入 `--primary-threshold=<实际默认层>`；工具默认 standard。完成 `report_review_checklist.csv` 后加 `--publish` 才能通过本地发布门禁，此选项不会部署到网上。

## 针对性再可视化(visualize_results.R)

`run_analysis.R` 跑完后，核心中间结果已存盘。`visualize_results.R` **读取保存结果做定制图，不重拟合 DESeq2/ORA/GSEA**，但会重新计算图内统计与布局。
脚本会切换到自身所在目录；使用项目本地副本，修改顶部 CONFIG 中的输入
绝对路径和独立输出目录。中央脚本的默认相对路径不适用于任意 run。
正式审阅优先用下文只读 RunReview notebook，避免改动已经归档的产物。

```bash
Rscript visualize_results.R
```

编辑文件顶部 CONFIG 块，三个独立 section（可分别开关）:

| section | 读取 | 用途 |
|---|---|---|
| `DO_KEY_GENES` | `DEG_results.Rdata`(vsd_mat) | 改 key genes → 单基因 bar/SEM + 热图 |
| `DO_GSEA_TERMS` | `gsea_results.rds`(gseaResult 对象） | 指定比较 + 条目（子串匹配或 top N)→ 逐条 gseaplot2 |
| `DO_ORA_THEME` | `2-GSEA/<threshold>/*.csv`(ORA 表） | 主题 dot-heatmap，可自定义主题字典 |

典型场景：跑完主流程后，想换一组 key genes 重画、或挑几条感兴趣的 GSEA 条目单独出出版级图、或按自己课题的主题（如"干扰素/抗原提呈/细胞毒性"）重排 ORA 主题图——改 CONFIG 重跑几秒即可，不用碰 notebook、不用重算。

`run_analysis.R` 已把 `gseaResult` 对象缓存到 `2-GSEA/gsea_results.rds`，所以单条 gseaplot2 不需要重跑 GSEA。旧项目若没有这个 `.rds`，`DO_GSEA_TERMS` 会提示并跳过；确需补算时新建 run 并记录父 run，不在旧目录重跑。

## RNAseq_lib 的定位

脚本按以下顺序寻找 `RNAseq_lib/`:
1. 环境变量 `RNASEQ_LIB_DIR`（最高优先级）
2. 调用目录或 runner 同目录下的 `RNAseq_lib`
3. runner 所在 git 仓库根目录的 `RNAseq_lib`（通过 rprojroot)
4. runner 上级目录的 `RNAseq_lib`

真实项目与 backend 分开保存，显式 `export RNASEQ_LIB_DIR=/path/to/pinned/backend/RNAseq_lib`，避免数据输出进入共享模板或随 backend 更新改变环境。

## 批量跑多个项目

每个项目分别使用 [BEST_PRACTICES](../../BEST_PRACTICES.md#1--reproducible-run)
中的新目录执行模式，并保存各自日志与退出状态。并行调度必须逐任务汇总
失败，不能仅凭后台任务已结束就宣布全部成功；随后检查必要产物和跳过的步骤。

## 什么时候回到 notebook?

- 想**换一组 gene set 重新做 GSVA** 并立刻看 heatmap/boxplot;
- 想对某个特定基因/通路做**一次性定制图**;
- 想逐 cell 检查中间对象（`vsd_mat`、`res_list`、`gsea_results`)。

`run_analysis.R` 已把 `vsd_matrix.csv`、`colData.csv`、`DEG_results.Rdata`、`step0-input.Rdata` 都存到 `1-DEG/`,notebook 或下游模板（WGCNA/TimeCourse）可直接读取这些中间结果，无需重跑。

General 提供两类 notebook：`RNAseq_General.ipynb` 是可独立完整执行的流程；
`RNAseq_General_RunReview.ipynb` 只读一个已完成的 `RUN_DIR`，把后续基因、
基因集和可视化调整写到独立 `REVIEW_OUTDIR`。真实项目常规生产以 CLI run
为统计事实来源，审阅 notebook 不修改该 run。
