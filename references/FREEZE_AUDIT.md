# 冻结审计与真实项目使用边界

更新日期：2026-08-31。维护版本：`v0.11.1`，基于 `v0.11.0`
（`9418686a9e110ae44dd00078ee3a63af4047351e`）完成代码、测试与文档修复。
本次没有升级全局依赖、创建 Git tag 或发布托管 release。原只读审计保留在
[历史发现](validation/2026-08-31/audit-findings.md)，不再代表当前未修复状态。

## 冻结判断

**可以作为有明确边界的稳定维护基线冻结；不能称为所有方法、平台和项目环境的完全验收。**
已关闭本次复现的实现缺陷，不需要继续扩张分析模块。真实项目仍须冻结自己的
输入、设计、后端 commit、环境和参考数据，并完成科学审阅。

| 层次 | 当前状态 | 含义 |
|---|---|---|
| 功能范围 | 冻结 | 保持六类生产 runner；只接受修复、兼容性、参考库适配和文档维护 |
| 本次缺陷修复 | 已验收 | 下表列出的代码和结果正确性回归已通过 |
| 本机生产主路径 | 已验证 | 完整测试、六模板 demo、额外人/鼠 TME CLI、General 报告与门禁通过 |
| 所有可选/网络/跨平台路径 | 不宣称完成 | IOBR 在线参考、GDC/GEO 下载、Quarto、旧 GSEA 后端等需要单独验证 |
| 长期环境复现 | 项目级责任 | 仓库没有已恢复验收的 lockfile/容器；版本清单不能代替可恢复环境 |

## 已关闭的维护项

| 原问题 | 修复 | 验收证据 |
|---|---|---|
| WGCNA hub 使用颜色查数字 eigengene，导出基因对基因相关 | 统一颜色命名；每基因只计算其所属模块 kME；按有效样本数算 module-trait P 值 | 多模块、grey、指定/空模块、缺失 trait 回归；完整 WGCNA demo 逐项复算 hub 值与唯一性 |
| 小鼠 native CIBERSORT 被错误转为人符号 | 保留原物种分支，人参考方法另建 HGNC ortholog 分支 | 人/鼠 native+ssGSEA 完整 CLI；547/547、511/511 参考基因重叠 |
| TPM 分母受到 DEG/biotype 筛选或提前合并符号影响 | 原始完整 feature universe 先算 TPM，后映射/去重；General 保留最终样本集 | TPM 列和、已知分母、样本排除及 native/human 映射回归 |
| 标为 ssGSEA 的结果实际调用 GSVA | 共享 `run_tme_ssgsea()` 调用 `ssgseaParam(normalize=TRUE)` | 与直接 GSVA ssGSEA 调用数值一致，输出 signature coverage |
| native 单样本、隐式对数转换、未置换伪 P 值 | 保留矩阵维度、显式输入尺度；perm=0 的 P 值写 NA；拒绝 VST/rlog | 单样本与低值线性/对数输入回归；runner/notebook 输入尺度拒绝测试 |
| native/IOBR 的不同参考被当等价比较 | 限 human、实际 reference/QN 一致，且输入最大值≥50，避免 IOBR 的低值自动反对数分支 | 不同物种、参考、QN、未知/低值输入均拒绝比较；不是独立生物学验证 |
| 仅主输入哈希/缺失校验和也被判重复 | schema-v2 增加 `run_inputs.csv` / `run_runtime.csv`；完整输入、分析、代码和运行环境签名才可参与候选分组 | 8 类辅助输入变化、NA/空/坏哈希、旧 schema、运行环境差异回归 |
| 复用目录覆盖配置与旧成果 | 六入口写文件前拒绝已有 native 输出；WGCNA 同时保护自定义子目录 | 六个真实 Rscript 子进程早停；目录文件及 MD5 保持不变 |
| limma `MIN_SAMPLE_FRAC` 未生效 | 原始 count 样本比例过滤与 edgeR `filterByExpr` 同时满足 | 不同比例/不平衡数据/按库大小过滤测试；0 保留旧 edgeR-only 行为 |
| GSEA 丢失不可计算通路 ID，且集合大小按原参考筛选 | 先按排名输入实际交集筛集合，避免上游 nominal-P 下标丢 ID；保留全部受检表及质量列 | 本机 GO、离线 KEGG 风格集合测试；BH 全家族复算、全 NA/部分 NA/零显著测试 |
| 非空 GSEA 表被当作显著/激活证据 | 自动图只用有效 BH 显著项；人工指定非显著曲线标明状态；报告排除无效行并显示数量 | 图与双报告模板回归、真实 HTML 渲染；不以 NES 方向断言激活/抑制 |
| 通路 bar 的默认 term 位置不一致 | 普通 ORA 保持条内；双向 ORA/GSEA 同步为零轴侧条内锚点，FDR 独立列 | 7 张单/双侧与长名称/短 bar 排版示例逐张渲染；不改变 bar 数值长度 |
| 专题报告在外部项目找不到共享模板 | 专题 runner 显式定位 backend report，保留 General 已有行为 | 仓库外、含空格工作目录下执行六入口实际报告分支回归 |

ssGSEA 与 GSVA 的参数类区别见 [GSVA 官方说明](https://bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.html)。
GSEA 兼容修复依据本机 `clusterProfiler 4.20.0` / `enrichit 0.2.0` 行为，并参考
[clusterProfiler 变更记录](https://github.com/YuLab-SMU/clusterProfiler/blob/devel/NEWS.md) 与
[enrichit 官方实现](https://raw.githubusercontent.com/YuLab-SMU/enrichit/master/R/gsea.R)。
未安装旧 DOSE/fgsea 后端，因此仅保留兼容分支，不声称该分支已实测。

## 最终验证

代码复制到隔离目录，添加仅用于资源定位的临时 Git 根；未覆盖仓库或真实项目
既有演示/分析输出。临时 Git 根没有发布 revision，测试源码由证据目录的
`SOURCE_SHA256SUMS.txt` 标识；最终采用版本应以本次维护提交为准。

| 检查 | 结果 |
|---|---|
| 完整 `tests/testthat.R` | **1,021 PASS，0 FAIL，0 SKIP，3 WARN**；告警来自 UpSetR/ggplot2 弃用接口 |
| notebook 与 R 语法 | 10 notebooks、17 共享 R 文件通过；最终仓库 85 个 R 文件（含验证脚本）解析通过 |
| 完整 smoke | General helper + 六个生产 runner 全部通过；KEGG 网络访问失败明确跳过，不能视作 KEGG 网络验收 |
| 人/鼠 native+ssGSEA CLI | 各 6 样本；22/25 细胞类，27/26 ssGSEA 签名；比例和=1，直接评分数值一致 |
| General 报告 | rmarkdown fallback 生成成功，草稿校验通过；未签署 `--publish` 按预期拒绝 coverage/scientific review |
| General 结果 QA | GO GSEA 1,513 受检/有效条目，0 无效、0 BH 显著；ID/实际交集/BH 全表复算一致 |
| PDF 与 HTML | 20 个 PDF 均可渲染、首页非空白；抽查 PCA、volcano、DEG heatmap、ORA 图；HTML 含内嵌预览 |
| 通用研究项目结构审计 | 0 error，15 missing-file warning；本仓库是软件后端，README/ROADMAP/CHANGELOG/本审计承担相应职责，未添加空研究项目文档 |

环境：macOS/aarch64，R 4.6.1、WGCNA 1.74、GSVA 2.6.3、babelgene 22.9。
[最终日志与校验脚本](validation/2026-08-31-maintenance/README.md) 和
[原版本证据](validation/2026-08-31/README.md) 分开保存。测试是工程与合成数据
验证，不是对任意真实队列的生物学签署，也不是所有图的投稿审阅。

## 历史项目迁移

1. 固定旧 run，记录旧 backend commit、配置和环境，保留旧结果。
2. 使用新 commit、新 `run_id`，填写 `PARENT_RUN_ID`、`RUN_CHANGE_NOTE`。
3. WGCNA hub 必须新算/重新导出正确 eigengene；native mouse、受 TPM 分母影响
   的 TME、旧 `ssGSEA_*` 和受影响 GSEA 必须重算。改名不能纠正算法。
4. limma 默认 0.5 现在实际生效，可能改变基因集合；若预先指定只用 edgeR 的
   组别感知过滤，明确设置 0。不要为了保留旧阳性结果事后调参。
5. 旧 manifest 不回写成 v2，不据旧 `duplicate_of` 删除或瘦身；重新登记后复核。
6. 更新所有引用旧结果的报告/图，保留 superseded 来源，再做科学审阅和归档。

## 真实项目的推荐用法

完整命令见 [BEST_PRACTICES](../BEST_PRACTICES.md)。不要另造 `md/template` 层：
当前文档职责是根目录入口、`templates/*` 配置/runner/专题说明、`references/*`
契约。项目只维护自己的 config、元数据、必要的定制脚本与审阅记录。

- **先设计**：记录生物学重复、物种、尺度、样本对应、排除规则、主比较、FDR/
  效应量阈值。DESeq2/voom 用未归一化 counts；TME 用 TPM 或 counts+匹配长度；
  WGCNA/聚类可用 VST/rlog/log2(TPM+1)。见 [DESeq2 官方输入与设计说明](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)。
- **固定后端和环境**：项目固定 commit；保存并实际验证环境恢复，以及 R、系统
  工具、注释/基因集/参考缓存。`sessionInfo.txt`/`run_runtime.csv` 不是恢复包；
  [renv 也有系统与外部资源边界](https://rstudio.github.io/renv/articles/renv.html)。
- **每版只保留一个生产执行来源**：在新 `analysis/runs/<run_id>/` 调中央 CLI，
  显式设置 `RNASEQ_LIB_DIR`。输入用绝对路径或相对 run 根，非相对 config；
  专题直接引用前 run 的输入，不把旧 `1-DEG/` 复制进新 run。完整 notebook 用于
  探索；已完成 General 用 RunReview 读取，衍生图另存。
- **验收所需范围**：核对实际模型、QC、有效结果、映射损失和图。必须的步骤失败
  或跳过就是未完成；关闭的不需要。TME mouse→human 仍是跨物种外推，本次
  synthetic mouse 有 90 基因未映射，成功执行不证明该人参考适合鼠组织。
- **报告到交付**：草稿→科学审阅→本地 publish gate→精选到 `results/`→归档。
  表、图、报告保留源 run/源表。`--publish` 不上传网站，也不能自动替人签署。

通用设计边界仍保留：General 的 `PAIR_ID` 自动使用 `~ pair_id + condition`，
不会保留额外 batch 项；复杂协变量/交互/嵌套模型需项目级适配。TimeCourse
按时间聚合可能合并 condition，提供 baseline 时间比较而非通用纵向混合模型。
运行目录保护不是并发锁；最终不可变归档也需项目自行实施。未知 KEGG/IOBR 等
参考将阻止候选重复分组，这是保护行为，不是分析失败。

下一步是按上述流程接入一个真实项目并验证其所需环境和方法，无需继续增加模板模块。
