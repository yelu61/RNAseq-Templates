# RNAseq_TME — command-line template

Non-interactive tumor-microenvironment (TME) deconvolution pipeline. It shares
the analysis intent and helper library of `notebooks/RNAseq_TME_Deconvolution_Template.ipynb`;
the production runner is authoritative: expression →
ESTIMATE / IOBR / native CIBERSORT → ssGSEA immune-signature scoring.

## Files

| File | Purpose |
|------|---------|
| `config.R` | **The only file you edit.** Paths, input mode, samples/groups, method switches. |
| `run_analysis.R` | Full pipeline: expression → TPM → ESTIMATE / IOBR / native CIBERSORT → ssGSEA → plots. |
| `visualize_results.R` | Cheap re-plotting from saved results (no deconvolution/ssGSEA recompute). |

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
next to the script, the repo root via `rprojroot`, or `../RNAseq_lib`). The
bundled CIBERSORT references are located via `rprojroot`
(`references/CIBERSORT/`) unless `CIBERSORT_SCRIPT` / `CIBERSORT_SIGNATURE` are
set explicitly in `config.R`.

## Input modes

- **`raw_counts` (default, recommended).** A featureCounts-like table with gene
  lengths (a length column, or start/end coordinates). Counts are converted to
  TPM with `extract_gene_lengths()` + `counts_to_tpm()` before deconvolution.
- **`expression`.** A pre-computed TPM or log2(TPM+1) matrix (`EXPR_UNIT`).
  VST/rlog is rejected — it cannot be converted back to a linear scale and is
  invalid for ESTIMATE/CIBERSORT/EPIC/xCell.

Mouse symbols are converted to human orthologs (TME methods use human
signatures); this needs network access for biomaRt/babelgene. Human-symbol input
can run offline when `RUN_IOBR=FALSE`, or after IOBR's method reference bundles
have been populated in the local cache.

## Output layout

Outputs use the numbered run layout (set `OUTDIR <- "."` in config.R to write
these at the project root; set a path to nest the run under a directory):

```
0-Config/analysis_config_used.R     frozen config snapshot (reproducibility)
4-TME/                              TPM_matrix.csv, ESTIMATE_scores.csv,
                                    IOBR_<method>.csv + IOBR_TME_combined.csv,
                                    CIBERSORT_native_*, ssGSEA_*, all PDFs,
                                    tme_results.Rdata (for visualize_results.R)
Analysis_summary.txt                text summary of the run
sessionInfo.txt                     package versions (reproducibility)
run_manifest.csv                    inputs/config/code signatures + lifecycle
```

## Optional methods and graceful degradation

ESTIMATE, IOBR, and native CIBERSORT are **optional**. Each is gated by a
`RUN_*` switch in `config.R`. A missing package or bundled native-CIBERSORT
resource disables that module with a clear message:

- `RUN_ESTIMATE` needs the `estimate` package.
- `RUN_IOBR` needs the `IOBR` package. On first use, recent IOBR releases may
  also fetch method reference bundles (`common_genes`, `TRef`, `xCell.data`,
  `lm22`). Populate that cache while online before an offline production run.
  If `RUN_IOBR=TRUE` and every requested method fails, the runner stops rather
  than silently reporting a complete deconvolution.
- `RUN_CIBERSORT` needs the bundled `references/CIBERSORT/CIBERSORT.R` script,
  a signature matrix (auto LM22 / `cibersort_mouse_22.csv`), and the
  `e1071`/`preprocessCore`/`future`/`furrr`/`purrr` packages.

ssGSEA (sections 8–9) needs `GSVA`; it too degrades gracefully if absent.

The bundled CI smoke test is deterministic and offline: it exercises native
ESTIMATE plus ssGSEA with `RUN_IOBR=FALSE`. IOBR should be validated separately
on a machine where its method data are already cached.

## Notes

- **Bug fix vs the notebook.** The notebook defines `group_df_for_plot` only in
  its visualization section, but the native-CIBERSORT section uses it earlier,
  so a top-to-bottom notebook run with `RUN_CIBERSORT = TRUE` fails. Here
  `group_df_for_plot` is built from the metadata right after loading, before any
  section that needs it.
- `SPECIES` selects the human vs mouse CIBERSORT signature and ortholog
  conversion automatically.
- `RUN_CIBERSORT_COMPARISON` (default TRUE) emits native-vs-IOBR CIBERSORT
  concordance tables and scatter/Bland–Altman PDFs when both variants run.
- The shared HTML report (`reports/analysis_report.qmd`) is oriented to the
  General (DESeq2) layout, so `GENERATE_HTML_REPORT` defaults to `FALSE`; the
  wiring is present if you add a TME-specific report template.
- `visualize_results.R` reloads `4-TME/tme_results.Rdata` and re-plots TME
  barplots/boxplots/heatmaps, the ESTIMATE boxplot, and the ssGSEA heatmap
  without re-running any deconvolution.
