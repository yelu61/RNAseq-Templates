# 冻结审计与实际项目使用边界

审计日期：2026-08-31。对象：`v0.11.0`，commit
`9418686a9e110ae44dd00078ee3a63af4047351e`。模式：只读实现审计、隔离执行验证、
有依据的文档维护；本次未修复分析源码、升级依赖、创建版本标签或发布 release。

## 结论

**可以继续维持 `v0.11.x` 功能冻结，但不能宣布全部验收完成或长期环境已经冻结。**
框架、生产入口、基础测试和报告机制已经成形；下一阶段应关闭维护缺陷，而非增加模块。
测试通过证明已覆盖路径能够执行，不等于每份科学产物正确。

| 冻结层次 | 当前判断 | 实际含义 |
|---|---|---|
| 功能范围与接口 | 可维持冻结 | 六类 runner 已有，停止扩张；允许修复、兼容性与文档维护 |
| 全路径科学产物验收 | 未完成 | WGCNA hub、小鼠 native CIBERSORT 有已复现缺陷 |
| 运行来源与留存验收 | 有条件 | manifest 可用，但辅助输入未完整校验，去重存在误判 |
| 长期环境复现 | 未冻结 | 没有仓库级 lockfile/容器；需项目自行保存可恢复环境和参考数据 |

## 本次验证记录

在 HEAD 的临时导出副本中测试，没有覆盖现有项目或仓库演示输出。副本加入了
临时 Git 根目录用于工具定位；分析源文件对应上述 commit，临时副本的 Git
revision 不用于发布溯源。现有测试依赖 Git 根目录，首次无 `.git` 的导出副本
出现资源定位失败；加入根目录后重跑，以下结果才是有效基线。

| 检查 | 本次结果 | 不能推导出的结论 |
|---|---|---|
| `Rscript tests/testthat.R` | 603 PASS，0 FAIL，0 SKIP，11 WARN | 不是零告警；不能覆盖未写断言的产物语义 |
| `Rscript tests/validate_notebooks.R` | 10 notebooks、17 共享 R 文件通过 | 只是语法/结构，不是逐 notebook 完整执行 |
| 逐一解析 Git 跟踪的 `.R` 文件 | 76 文件通过 | 不保证运行时依赖或计算结果 |
| `Rscript examples/run_demo_smoke_test.R` | General helper + 六个生产 runner 演示通过 | 网络 KEGG 步骤可跳过；IOBR/native CIBERSORT 未在 TME demo 开启 |
| General 完整 demo 的独立报告渲染 | rmarkdown fallback 渲染及草稿验证成功 | 不是 Quarto 分支验收，也不是科学签署 |
| 未签署草稿执行 `--publish` | 按预期拒绝：coverage、scientific review 未完成 | 没有替用户签署，不宣称 publish-ready |
| 报告结构校验 | HTML/内嵌预览、DEG 数量复算及 30 个 PDF 非空检查通过 | 非逐图人工视觉或生物学审阅 |

环境：R 4.6.1，macOS/aarch64；WGCNA 1.74、babelgene 22.9。
11 个单元测试告警包括 UpSetR/ggplot2 弃用提示和 GSEA 无法计算部分通路、
NA 标识等提示。请保留并审阅，不把它们等同于全部科学验收通过。
完整日志、三个最小缺陷复现和 GSEA 输出检查见 [验证证据目录](README.md)。
未查询当前远端 CI 状态；这里报告的是本机本次执行，不是历史徽章状态。

## 未关闭的维护项

**CIBERSORT 范围澄清：** `references/CIBERSORT/CIBERSORT.R` 是不依赖 IOBR 的
本地实现；`cibersort_mouse_22.csv` 和 `LM22.txt` 分别提供小鼠与人参考。
下表的问题限定于 TME runner 的输入接线，不是本地算法缺失或小鼠参考不可用。
直接调用 native 时应保留 MGI 输入并显式选择小鼠 CSV；只有需要人类参考的
其他方法才使用相应的人源输入。参见 [本地 CIBERSORT 说明](../../CIBERSORT/README.md)。
人/鼠输入的推荐转换阶段、各方法分支及两个入口的现有差异，见
[TME 转换与分流契约](../../../templates/TME/README.md#species-and-symbol-routing-required-contract)。

| 优先级 / 问题 | 证据及影响 | 修复前使用限制 | 验收条件 |
|---|---|---|---|
| P1：WGCNA hub 导出使用错误 eigengene | [runner](../../../templates/WGCNA/run_analysis.R) 使用数字标签的 ME1/ME2，却按颜色找 MEturquoise；NULL 使相关运算变成基因对基因。最小例 60 基因导出 3,600 行、3,540 重复；完整 demo 338,452 行仅 2,986 唯一基因，仍通过 smoke | 不解释 `Hub_genes_*`，不能声称该路径已验收；网络和 module-trait 结果仍需独立 QA | 统一 module/eigengene 映射；每基因一行，kME 与其真实模块 eigengene 的直接计算一致；多模块、grey、定向模块选择均覆盖 |
| P1：小鼠 native CIBERSORT 物种错配 | [TME runner](../../../templates/TME/run_analysis.R) 先将表达转 HGNC，再按 `SPECIES=mouse` 选择 MGI signature；511 行参考构造的输入转换为 404 行，精确重叠仅 3 行 | 小鼠专题 run 关闭 `RUN_CIBERSORT`；不将其他方法的成功视作此分支成功 | 明确每个方法表达/参考的物种与 ID 空间；增加 human/mouse 端到端重叠和输出检查，不能只改大小写 |
| P1：候选重复被误称“完全重复” | [run_utils.R](../../../RNAseq_lib/run_utils.R) 只校验一个主输入；同路径 metadata 改内容仍被分组；缺失 input MD5 的在线 run 也会被分组，两项均已复现 | `duplicate_of` 只用于提醒复核，不能据此删除/瘦身；项目另存所有输入与参考文件校验和 | metadata/traits/clinical/TPM/annotation/cache 内容变化须影响签名；NA/空校验和不得证明相同；环境与参考版本相等另行核对 |
| P2：运行目录没有防覆盖保护 | runner 会复用目录、覆盖 config/manifest，并可能留下上次可选模块成果 | 每次运行必须新建目录；保留失败运行，修订后换 run ID | 用 [安全执行示例](../../../BEST_PRACTICES.md#1--reproducible-run) 拒绝已存在目录；若增加程序防护，应测试旧目录拒绝且内容不变 |
| P2：limma `MIN_SAMPLE_FRAC` 未生效 | [helper](../../../RNAseq_lib/limma_voom_utils.R) 接收此参数但 `filterByExpr()` 不使用它 | 不把该参数写成已实施过滤条件 | 明确预期过滤语义后实现或废弃，增加能区分参数取值的输入检查 |
| P1：TME 免疫评分算法名称与调用不符 | 追查 symbol 分支时发现，专题 TME 和 General 的“ssGSEA”段均调用 `gsvaParam()`，实际选择 GSVA；[官方文档](https://bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.html) 明确区分它与 `ssgseaParam()` | 不把现有 `ssGSEA_*` 文件直接写成已验证的 ssGSEA 方法结果 | 确认预期算法，纠正参数类/命名，增加与直接 ssGSEA 调用的一致性回归；历史结果以新 run 重算或明确更正方法 |
| P2：鼠 native 与人参考 IOBR 的分类不等价 | 当前自动比较未限定参考物种；鼠 CSV 25 类，人 LM22 22 类；清理列名不是细胞类型本体映射 | 小鼠项目关闭自动 `RUN_CIBERSORT_COMPARISON`，不能直接把两套比例当作同量纲验证 | 比较前验证物种、参考及分类一致；不同参考仅允许预定义且有记录的探索性大类比较 |
| P2：当前依赖下 GSEA 表存在不可用行，排名图不代表显著 | 本次 General CLI 的 GO 表 2,337 行，2,028 行全 NA；剩余 309 行中 BH<0.05 为 0。本地 `clusterProfiler 4.20.0` / `enrichit 0.2.0` 的 P 值逻辑筛选未排除 NA，使不可计算条目丢失原 ID；runner 仍按非空表绘图 | 逐项核查有效 ID、有限统计量、FDR 与不可计算项；非显著 top term 只能描述为排名展示，不称通路激活/抑制证据 | 验证并固定兼容依赖，保留不可计算项及原因；表/图/报告明确有效和显著条目数，覆盖全 NA、部分 NA、零显著三种情形 |

上述 P1 是“全路径验收完成”的阻碍，但不要求继续增加功能。修复涉及真实统计
结果的旧项目，应先识别是否走过受影响路径，再以新 run 复算；不覆盖原报告。

## 使用范围与环境限制

- **General / Limma_Voom**：可作为预先审查设计后的常规组间差异分析入口。
  General 的 `PAIR_ID` 自动设为 `~ pair_id + condition`，会覆盖手写公式；
  同时设置 `BATCH_VECTOR` 不等于同时调整批次。任意协变量、复杂配对、交互
  和嵌套设计不能仅凭配置名称就宣称已支持。图上的简单两两检验也不能替代模型。
- **TimeCourse**：Mfuzz 按时间聚合，可能合并不同 condition；DEG 是相对
  baseline 的时间效应，配对使用固定 subject 项。它不是通用纵向混合模型或
  treatment×time 交互检验器。
- **TME / 公共数据**：原始 counts 加匹配的 feature/exon 长度，或已核实 TPM；
  不把 VST/rlog 当 TPM。方法参考缓存、物种映射、数据库版本及网络失败需单独
  记录。本地 cohort 文件并不意味着 KEGG 等步骤完全离线。
- **报告**：General 有标准 HTML contract；专题模板默认关闭此报告。CI 的
  demos 均关闭 HTML；本次另验了 General 的 rmarkdown 路径，未验所有可选
  模块、在线下载、操作系统或逐图科学正确性。
- **GSEA**：本次表中的 2,028 个全空行不是 2,028 个有效通路。报告会滤除
  全空行后展示 top 条目，但并不要求这些条目达到预设 FDR；图上 NES 方向也
  不能替代显著性判断。删除空行只能改善产物格式，不能把本次无显著结果变为
  科学阳性结论，亦不能恢复丢失的通路标识。
  309 是经过 nominal-P 筛选的有限导出行数，不是全部受检通路总数；验收还需
  检查基因集与输入的实际重叠，不能只看配置中的最小集合大小。
- **环境**：CI 使用浮动的 `ubuntu-latest`、R `release` 与未固定到 revision
  的依赖来源。`sessionInfo.txt` 是记录，不是恢复包。`dorothea`/`viper` 尚未
  被默认安装器列入，开启 TF 需额外预检。lockfile 也不保存 R 本体、系统工具
  或在线参考数据；参见 [renv 官方边界](https://rstudio.github.io/renv/articles/renv.html)。

## 真实项目的推荐操作

遵循 [BEST_PRACTICES](../../../BEST_PRACTICES.md)，不再新增一套重复模板目录：

1. 先记录问题、统计单位、样本对应、输入尺度、主比较和主阈值，再配置分析。
2. 固定 backend commit 与项目环境。项目仅维护自己的 config、元数据、基因集
   和必要的定制脚本；共用模板的修改必须走维护验证。
3. 一个确定的分析版本只保留一个生产执行来源：CLI 写入新的
   `analysis/runs/<run_id>/`。完整 notebook 仅供探索；已完成 CLI 的跟进使用
   RunReview 只读源 run，衍生图写到独立位置。
4. 必须完成的步骤不能意外跳过；查看 QC、实际模型、有限统计量、基因集映射、
   输出表和图。必要的敏感性分析另建 run，说明与主 run 的关系。
5. 草稿报告 → 科学审阅 → 本地 publish gate → 精选到 `results/`；所有交付
   文件都保留源 run、源表、图参数和审阅记录，最终归档后再视为不可变。

## 文档维护与下一步

本仓库没有 `md/template` 目录。文档职责分别是：根目录 README 负责入口，
BEST_PRACTICES 负责项目生命周期，GETTING_STARTED 负责教学，`templates/*/README.md`
负责各路径限制，`references/` 负责参数/产物/审计契约。此次修正了这些文档中
已验证的陈旧表述；CHANGELOG 保留历史发布记录，ROADMAP 增加未关闭维护项。

通用项目管理审计显示 0 error、9 项研究项目文件缺失 warning；本仓库是软件
后端，README/ROADMAP/CHANGELOG/本审计承担对应职责，因此未新增一套无实际
样本数据的 PROJECT_BRIEF、DATA_DICTIONARY 等研究项目空壳文件。

下一最小任务：先修复 WGCNA hub 导出并增加正确性回归，再处理 mouse native
CIBERSORT 与完整输入签名；完成相关验证后更新本审计，不以本次文档修订冒充
bug 已修复。可用以下交接请求继续：

> 基于 references/FREEZE_AUDIT.md，先修复 WGCNA hub 的 eigengene 映射与
> 导出正确性，为现有错误添加回归测试，再运行完整测试和六模板 smoke。
> 保持 v0.11.x 功能冻结，不增加新模块，不覆盖真实项目历史 run。
