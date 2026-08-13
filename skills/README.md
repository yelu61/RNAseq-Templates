# Skills

本仓库**不再内置** agent skill。统一的 `bulk-rnaseq-analysis` skill 已迁移到中央集合维护，本仓库仅作为其后端之一被路由调用。

## 规范位置（canonical）

- GitHub：<https://github.com/yelu61/agent-ready-research-skills> → `skills/bulk-rnaseq-analysis`
- 本地克隆：`~/Projects/agent-ready-research-skills/skills/bulk-rnaseq-analysis`
- 用户级发现链接：`~/.agents/skills/bulk-rnaseq-analysis`（→ 中央集合），并经 `~/.claude/skills/bulk-rnaseq-analysis` 暴露给 Claude Code

## 使用

skill 在中央集合安装/链接后即可全局调用，无需在本仓库内放置任何文件：

```
Use $bulk-rnaseq-analysis to inspect these counts and metadata ...
```

它会按显式后端契约把"本地/GEO 矩阵类"工作路由到本仓库（RNAseq-Templates），把"TCGA/TARGET/GTEx、MAF/TMB、CNV、甲基化、预后、药敏"等工作路由到 TCGA toolkit；它不复制、不合并任一端的实现。

## 为什么仓库里不再有 skills/

此前 `skills/bulk-rnaseq-analysis` 曾短暂内置于此，后迁至中央集合统一演进。仓库内 `.agents/skills/` 与 `.claude/skills/` 下指向旧路径的符号链接已随迁移清理（目标已删除，原为悬空链接）。如需改动 skill，请到中央集合提交，而非本仓库。
