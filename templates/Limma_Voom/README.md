# RNAseq_Limma_Voom — command-line template

Non-interactive limma-voom differential-expression pipeline. It shares the
analysis intent and helper library of `notebooks/RNAseq_limma_voom_Template.ipynb`;
the production runner is authoritative and may harden ordering or output layout.

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
run_manifest.csv                    inputs/config/code signatures + lifecycle
```

## When to use limma-voom vs the General (DESeq2) template

Both answer "which genes differ between groups" from raw counts. Choose the
method in the analysis plan, based on the study design and required workflow;
sample count alone does not determine which method to use. This runner exposes
group contrasts and an optional batch covariate, not every model supported by
the underlying limma package. General additionally provides the multi-threshold
DEG grid, GSVA, optional TME/TF and the standard HTML-report workflow.

Genes must pass both edgeR's group-aware `filterByExpr(min.count = MIN_COUNT)`
and the explicit raw-count rule: counts ≥ `MIN_COUNT` in at least
`ceiling(MIN_SAMPLE_FRAC * total_samples)` samples. `MIN_SAMPLE_FRAC` is in
`[0, 1]`; setting it to `0` keeps edgeR filtering alone. It is not edgeR's
`min.prop`, which has a different role in large-group filtering; see the
[edgeR reference manual](https://bioconductor.org/packages/release/bioc/manuals/edgeR/man/edgeR.pdf).

**This fixes a previously ignored parameter and can change results.** The
default `0.5` now enforces the configured total-sample fraction. In an
unbalanced design, that can remove genes expressed only in a small group;
set `0` when the study calls for edgeR's group-aware filtering alone, or to
reproduce the previous filtering behavior. Record the choice before analysis
and re-run downstream fitting and enrichment when filtering changes.

The runner rejects existing native output artifacts before writing a new
configuration snapshot. Completed runs also contain `run_inputs.csv` and
`run_runtime.csv`; unknown references make a run ineligible for duplicate
grouping. These files record provenance, not a restorable environment.

## Notes

- `SPECIES` selects the OrgDb (`org.Hs.eg.db` / `org.Mm.eg.db`) and KEGG organism
  code automatically — no hard-coded annotation package.
- `BATCH_VECTOR` (optional) is added as a covariate to the limma design.
- Single-term GSEA figures that are degenerate (too few genes) are skipped with a
  message rather than aborting the run.
- The shared HTML report (`reports/analysis_report.qmd`) is oriented to the
  General layout, so `GENERATE_HTML_REPORT` defaults to `FALSE` here; the wiring
  is present if you add a limma-specific report template.
