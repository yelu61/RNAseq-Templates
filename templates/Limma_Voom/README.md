# RNAseq_Limma_Voom — command-line template

Non-interactive limma-voom differential-expression pipeline. This is the
headless counterpart of `notebooks/RNAseq_limma_voom_Template.ipynb` — same
steps, driven by a config file instead of a parameter cell.

## Files

| File | Purpose |
|------|---------|
| `config.R` | **The only file you edit.** Paths, samples, groups, comparisons, thresholds. |
| `run_analysis.R` | Full pipeline: counts → voom → contrasts → volcano/heatmap → ORA + GSEA → theme maps. |
| `visualize_results.R` | Cheap re-plotting from saved results (no voom/GSEA recompute). |

## Usage

```bash
# 1. Copy this folder into your project, then edit config.R.
# 2. Run the standard analysis:
Rscript run_analysis.R                    # uses ./config.R
Rscript run_analysis.R path/to/config.R   # uses a different config

# 3. (Optional) iterate on figures without re-running the pipeline:
Rscript visualize_results.R
```

The runner preserves the invocation directory. Relative inputs and `OUTDIR`
resolve from that run root, while a relative config argument is resolved before
loading. For production, invoke it from a fresh `analysis/runs/<run_id>/`.

`RNAseq_lib` is located automatically (env `RNASEQ_LIB_DIR`, a `./RNAseq_lib`
next to the script, the repo root via `rprojroot`, or `../RNAseq_lib`).

## Output layout

Outputs use the numbered run layout (set `OUTDIR <- "."` in config.R to write
these at the project root; set a path to nest the run under a directory):

```
0-Config/analysis_config_used.R     frozen config snapshot (reproducibility)
1-DEG/                              per-comparison DEG tables + limma_voom_results.Rdata
2-GSEA/                             GO/KEGG ORA + GSEA result tables
3-Visualization/                    voom trend, volcano, heatmap, ORA/GSEA suites,
                                    ThemeEnrichment/ (theme dot-heatmaps, single-term gseaplot2)
Analysis_summary.txt                text summary of the run
sessionInfo.txt                     package versions (reproducibility)
```

## When to use limma-voom vs the General (DESeq2) template

Both answer "which genes differ between groups". Prefer **limma-voom** for
larger sample sizes, when you want limma's empirical-Bayes moderation and
flexible contrast/design matrices, or to stay consistent with a microarray-era
workflow. Prefer **General (DESeq2)** for small-n studies where count-based
dispersion shrinkage is more conservative, and when you want the multi-threshold
DEG grid, GSVA, TME and TF options that the General runner ships with.

## Notes

- `SPECIES` selects the OrgDb (`org.Hs.eg.db` / `org.Mm.eg.db`) and KEGG organism
  code automatically — no hard-coded annotation package.
- `BATCH_VECTOR` (optional) is added as a covariate to the limma design.
- Single-term GSEA figures that are degenerate (too few genes) are skipped with a
  message rather than aborting the run.
- The shared HTML report (`reports/analysis_report.qmd`) is oriented to the
  General layout, so `GENERATE_HTML_REPORT` defaults to `FALSE` here; the wiring
  is present if you add a limma-specific report template.
