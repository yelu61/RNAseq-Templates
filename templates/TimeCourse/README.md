# RNAseq_TimeCourse — command-line template

Non-interactive time-course (Mfuzz) pipeline. It shares the analysis intent and
helper library of `notebooks/RNAseq_TimeCourse_Template.ipynb`; the production
runner is authoritative:

1. Aggregate normalized expression to per-time-point means.
2. Mfuzz soft clustering (+ core-gene heatmap, trend panels, per-cluster GO ORA,
   and a publication-grade theme dot-heatmap).
3. Optional time-point-vs-baseline DESeq2 on raw counts (per-timepoint DEG tables
   + volcano + summary barplot).

## Files

| File | Purpose |
|------|---------|
| `config.R` | **The only file you edit.** Paths, time/group columns, Mfuzz + DEG settings. |
| `run_analysis.R` | Full pipeline: expression → time-point means → Mfuzz clustering + ORA → timepoint DEG. |
| `visualize_results.R` | Cheap re-plotting from saved results (no Mfuzz/DESeq2 recompute). |

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

## Two inputs, two jobs

This template reads **two** inputs, matching the notebook:

- **Normalized expression** (`EXPR_FILE`, e.g. a VST matrix from the General
  template) + `META_FILE` (with `sample`, `time`, optional `condition`) — used
  for aggregation and Mfuzz clustering. `validate_expression_contract` checks it
  is normalized (not raw counts).
- **Raw integer counts** (`RAW_COUNTS_FILE`) + `COUNT_META_FILE` — used only for
  the time-point-vs-baseline DESeq2 (`RUN_TIMEPOINT_DEG`). If
  `RUN_TIMEPOINT_DEG` is `TRUE` and the raw-counts file is missing, the run
  **hard-stops** with a clear message (set `RUN_TIMEPOINT_DEG <- FALSE` to skip).

## Output layout

Outputs use the numbered run layout (set `OUTDIR <- "."` in config.R to write
these at the project root; set a path to nest the run under a directory):

```
0-Config/analysis_config_used.R     frozen config snapshot (reproducibility)
1-DEG_Timepoint/                    per-timepoint DEG tables + Timepoint_DEG_summary.csv
3-Visualization/                    volcano per timepoint, Timepoint_DEG_summary.pdf,
                                    ThemeEnrichment/ (Mfuzz cluster theme dot-heatmap)
5-TimeCourse/                       mean_expression_by_time.csv, mfuzz_clusters.csv,
                                    mfuzz_trends.pdf, mfuzz_core_heatmap.pdf,
                                    GO_ORA_Cluster_*.csv, timecourse_results.Rdata
Analysis_summary.txt                text summary of the run
sessionInfo.txt                     package versions (reproducibility)
run_manifest.csv                    inputs/config/code signatures + lifecycle
```

Mfuzz-specific artefacts are grouped under `5-TimeCourse/` so the standard
numbered areas (`1-`, `3-`) stay comparable with the other templates; the saved
`timecourse_results.Rdata` (Mfuzz object/clusters, mean-by-time matrix, timepoint
DEG results, config values) lets `visualize_results.R` re-plot without recompute.
The trend figure uses low-opacity within-cluster trajectories plus a weighted
cluster centroid, avoiding Mfuzz's saturated rainbow default while retaining
the underlying temporal heterogeneity.

## Notes

- `SPECIES` selects the OrgDb (`org.Hs.eg.db` / `org.Mm.eg.db`) for the Mfuzz
  cluster ORA — fixing a notebook bug that hard-coded `org.Hs.eg.db`.
- `MFUZZ_SEED` is applied via `set.seed()` so Mfuzz clustering is reproducible.
- Mfuzz is gated behind `RUN_MFUZZ` **and** package availability: if the `Mfuzz`
  package is not installed the runner downgrades with a `message()` and still
  produces the time-point DEG, rather than aborting (mirrors the `RUN_TME`
  pattern in `templates/General/`).
- A single-level `condition` column is dropped from the DESeq2 design (a
  constant covariate breaks the model); multi-level conditions are included.
- `SUBJECT_COL` (optional) enables a paired/repeated-measures design
  (`~ subject + time`).
- The shared HTML report (`reports/analysis_report.qmd`) is oriented to the
  General layout, so `GENERATE_HTML_REPORT` defaults to `FALSE` here; the wiring
  is present if you add a time-course-specific report template.

## When to use this template

Use **TimeCourse** when your samples span an ordered time gradient and the
question is "which genes share a temporal expression pattern" (Mfuzz) and/or
"which genes change at each time point relative to baseline" (DESeq2). For a
simple two-group comparison use **Limma_Voom** or **General (DESeq2)** instead.
