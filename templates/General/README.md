# General 模板 — 命令行运行（run_analysis.R + config.R)

`notebooks/RNAseq_General.ipynb` 适合**交互式探索**（改 gene set 做 GSVA、临时加可视化、逐 cell 检查）。但在实际项目里，标准流程更适合**一键、可复现、可批量**地跑完——这就是 `run_analysis.R` + `config.R` 的用途。两者共享分析意图与 helper；生产 runner 是最终统计事实来源，产物不要求逐文件一致。

## 用法

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
   Rscript run_analysis.R                 # 用当前目录的 config.R
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
多个 run 后执行 `Rscript tools/build_run_registry.R analysis/runs`，统一生成
`RUN_REGISTRY.csv` 并识别完全重复 run；工具只报告，不自动清理文件。

## 可选分析开关(config.R)

| 开关 | 作用 | 依赖 |
|---|---|---|
| `RUN_TF_ANALYSIS` | DoRothEA/VIPER 转录因子活性 | dorothea, viper, limma |
| `RUN_TME` | 肿瘤微环境去卷积（ESTIMATE / IOBR / ssGSEA)，输出到 `4-TME/` | IOBR, estimate；小鼠需 biomaRt（联网） |
| `RUN_COMPARECLUSTER` | 多组 compareCluster 富集比较 | ≥3 组才有意义 |
| `EXPORT_EXCEL` | 导出多阈值 DEG 到 `1-DEG/DEG_results.xlsx` | openxlsx |
| `GENERATE_HTML_REPORT` | 渲染统一 HTML 报告 | quarto CLI 或 rmarkdown |

`THRESHOLD_GRID` 可为某一行添加 `p_column = "pvalue"`，作为低 DEG 情形下补充 ORA 的**探索层**；默认层的图、交集和报告会自动使用该层实际解析出的 P 值列。名义 P 层不能作为主结论阈值。

**TME 特别注意**：去卷积用 TPM，而标准流程用 VST，两者不可逆转换。所以 `RUN_TME=TRUE` 时脚本会从原始 counts + 基因长度重新算 TPM——你的注释文件需要基因长度信息：
- 有长度列 → 设 `TME_GENE_LENGTH_COLUMN <- "gene_length"`;
- 否则用起止坐标 → 设 `TME_GENE_START_COL` / `TME_GENE_END_COL`(featureCounts 注释通常自带）。
- 小鼠数据会自动经 biomaRt 转人源直系同源（需联网）；人源数据离线即可。
- 小鼠项目建议设置 `TME_ORTHOLOG_CACHE` 到项目的 `data/processed/`：首次成功映射后，重复运行不再依赖网络；已有缓存但新增 symbol 联网失败时会保留已缓存的可映射基因并给出 warning。
- General runner 默认 IOBR 方法为 `estimate/cibersort/epic/mcpcounter`；`xcell` 需要在线 biomaRt，只有在网络和参考库已验证时才显式加入。
- 子开关 `RUN_TME_ESTIMATE` / `RUN_TME_IOBR` / `RUN_TME_SSGSEA` 可分别控制三个方法。

HTML 报告的项目说明不要写进共享模板或已生成 HTML。首次运行时设 `SCAFFOLD_REPORT_INTERPRETATION <- TRUE`，再编辑 `report_interpretation/<section>.md`；这些文件会在重渲染时保留并嵌入相应章节。运行 `Rscript tools/render_report.R <run_dir> --scaffold` 可独立重渲染；完成 `report_review_checklist.csv` 后加 `--publish` 才能通过发布门禁。

## 针对性再可视化(visualize_results.R)

`run_analysis.R` 跑完后，所有中间结果都已存盘。`visualize_results.R` **直接读这些结果做定制图，不重跑任何计算**:

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

`run_analysis.R` 已把 `gseaResult` 对象缓存到 `2-GSEA/gsea_results.rds`，所以单条 gseaplot2 不需要重跑 GSEA（借鉴 202604DYY 的缓存思路）。旧项目若没有这个 `.rds`,`DO_GSEA_TERMS` 会提示并跳过，重跑一次 `run_analysis.R` 即可生成。

## RNAseq_lib 的定位

脚本按以下顺序寻找 `RNAseq_lib/`:
1. 环境变量 `RNASEQ_LIB_DIR`（最高优先级）
2. 当前目录下的 `./RNAseq_lib`
3. 所在 git 仓库根目录的 `RNAseq_lib`（通过 rprojroot)
4. 上级目录的 `../RNAseq_lib`

最省事的做法是把项目放在仓库克隆里（脚本会自动找到仓库的 `RNAseq_lib`)，或显式 `export RNASEQ_LIB_DIR=/path/to/RNAseq_lib`。

## 批量跑多个项目

```bash
for proj in projA projB projC; do
  (cd "$proj" && Rscript run_analysis.R > run.log 2>&1) &
done
wait
echo "all done"
```

## 什么时候回到 notebook?

- 想**换一组 gene set 重新做 GSVA** 并立刻看 heatmap/boxplot;
- 想对某个特定基因/通路做**一次性定制图**;
- 想逐 cell 检查中间对象（`vsd_mat`、`res_list`、`gsea_results`)。

`run_analysis.R` 已把 `vsd_matrix.csv`、`colData.csv`、`DEG_results.Rdata`、`step0-input.Rdata` 都存到 `1-DEG/`,notebook 或下游模板（WGCNA/TimeCourse）可直接读取这些中间结果，无需重跑。

General 提供两类 notebook：`RNAseq_General.ipynb` 是可独立完整执行的流程；
`RNAseq_General_RunReview.ipynb` 只读一个已完成的 `RUN_DIR`，把后续基因、
基因集和可视化调整写到独立 `REVIEW_OUTDIR`。真实项目常规生产以 CLI run
为统计事实来源，审阅 notebook 不修改该 run。
