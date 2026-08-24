# RNAseq_WGCNA — command-line template

Non-interactive WGCNA co-expression pipeline. It shares the analysis intent and
helper library of `notebooks/RNAseq_WGCNA_Template.ipynb`; the production
runner is authoritative. It builds co-expression modules from a
normalized expression matrix, correlates module eigengenes with sample traits,
and exports per-module hub genes.

There is **no `wgcna_utils` module** — the runner calls the `WGCNA` package
directly and sources only `plot_utils.R`, `data_utils.R` and `report_utils.R`
from `RNAseq_lib`.

## Files

| File | Purpose |
|------|---------|
| `config.R` | **The only file you edit.** Expression/trait paths, MAD filter, network and merge parameters. |
| `run_analysis.R` | Full pipeline: load → MAD filter / sample QC → soft-threshold → blockwise modules → module-trait correlation → hub-gene export. |
| `visualize_results.R` | Cheap re-plotting from saved results (no `blockwiseModules` recompute). |

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

## Input expectations

- `EXPR_FILE`: genes × samples **normalized** expression — VST / rlog /
  log2(TPM+1). Do **not** use raw counts. The runner validates this
  (`validate_expression_contract(expected = "vst")`).
- `TRAIT_FILE`: sample metadata containing `SAMPLE_COLUMN`. Sample-ID mirror
  columns are excluded; varying numeric traits are retained and categorical
  traits are encoded as labelled non-reference indicators (for example,
  `condition: Treatment`). All-NA and invariant traits are dropped. If no
  usable trait remains, the run stops with a clear message.

## Output layout

WGCNA tables/figures are written under `OUTDIR` (default `5-WGCNA`); the config
snapshot, run summary and session info are written next to it at the run root:

```
0-Config/analysis_config_used.R   frozen config snapshot (reproducibility)
5-WGCNA/
  Sample_clustering.pdf           sample dendrogram (QC)
  Soft_threshold_selection.pdf    scale-free fit + mean connectivity
  Module_dendrogram.pdf           gene dendrogram + module colors
  Module_trait_heatmap.pdf        module eigengene vs trait correlation
  WGCNA_gene_modules.csv          gene -> module color assignment
  WGCNA_network.rds               net, MEs, moduleColors, datExpr, traits
  Module_trait_correlation.csv    module-trait Pearson correlations
  Module_trait_pvalue.csv         corresponding Student p-values
  Hub_genes_<module>.csv          per-module hub genes ranked by |kME|
  Hub_genes_all_modules.csv       combined hub-gene table
  WGCNA_results.Rdata             targeted save for visualize_results.R
Analysis_summary.txt              text summary of the run
sessionInfo.txt                   package versions (reproducibility)
run_manifest.csv                  inputs/config/code signatures + lifecycle
```

## Key parameters

- `MIN_MAD_QUANTILE`: keep genes with MAD ≥ this quantile (top variable genes).
- `NETWORK_TYPE`: `"signed"` (recommended), `"unsigned"`, or `"signed hybrid"`.
- `SOFT_POWER`: soft-threshold override. `NULL` = auto-pick via
  `pickSoftThreshold` (falls back to the power with max scale-free R² when no
  automatic estimate is returned). Set an integer to force a power.
- `MIN_MODULE_SIZE`, `MERGE_CUT_HEIGHT`: module detection / merging.
- `TARGET_MODULES`: restrict hub-gene export to these module colors (`NULL` = all).

## Notes

- Soft-threshold selection is a judgment checkpoint: the runner uses
  `sft$powerEstimate`, and if that is `NA` it falls back to the power with the
  maximum signed scale-free fit and messages the choice. Review
  `Soft_threshold_selection.pdf` before trusting the modules.
- The shared HTML report (`reports/analysis_report.qmd`) is oriented to the
  General (DESeq2) layout, so `GENERATE_HTML_REPORT` defaults to `FALSE` here;
  the wiring is present if you add a WGCNA-specific report template.
