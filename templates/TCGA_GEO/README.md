# RNAseq_TCGA_GEO — command-line template

Non-interactive TCGA/GEO public-data mining pipeline. This is the headless
counterpart of `notebooks/RNAseq_TCGA_GEO_Template.ipynb` — same steps, driven
by a config file instead of a parameter cell.

It acquires a Tumor/Normal cohort three ways (GDC download, GEO SeriesMatrix, or
local files), cleans the clinical table, runs DESeq2 Tumor-vs-Normal, then does
single-gene expression, survival (KM by median/quantile, clinical KM, uni- and
multivariate Cox), and ORA + GSEA with publication-grade theme maps.

## Files

| File | Purpose |
|------|---------|
| `config.R` | **The only file you edit.** Data source, paths, thresholds, survival genes. |
| `run_analysis.R` | Full pipeline: acquire → clinical cleanup → DESeq2 → expression + survival → ORA + GSEA. |
| `visualize_results.R` | Cheap re-plotting from saved results (no DESeq2/GSEA/survival recompute). |

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

## Data source (pick exactly one in `config.R`)

| Mode | Flag(s) | Network | Notes |
|------|---------|---------|-------|
| GDC download | `DOWNLOAD_FROM_GDC <- TRUE` | required | TCGAbiolinks `GDCquery`/`GDCdownload`/`GDCprepare`; assay auto-detected by name (`count|unstrand`, `tpm`). |
| GEO download | `DOWNLOAD_FROM_GEO <- TRUE` | required | GEOquery SeriesMatrix via `GEO_ACCESSION`. |
| **Local files** | both `FALSE` | **none** | Reads `LOCAL_COUNTS_FILE` / `LOCAL_TPM_FILE` / `LOCAL_CLINICAL_FILE`. Fully offline — this is what the demo uses. |

The GDC/GEO branches are gated behind `requireNamespace("TCGAbiolinks")` /
`requireNamespace("GEOquery")` and stop with a clear message if the package (or
network) is unavailable. DESeq2 requires raw integer counts; a GEO SeriesMatrix
that is normalized/log expression is rejected with a pointer to a limma workflow.

## Output layout

Outputs use the numbered run layout (set `OUTDIR <- "."` in config.R to write
these at the project root; set a path to nest the run under a directory). The
notebook wrote files flat; this runner reconciles them into numbered dirs:

```
0-Config/analysis_config_used.R      frozen config snapshot (reproducibility)
1-DEG/DESeq2_Tumor_vs_Normal.csv     Tumor-vs-Normal DEG table
1-DEG/TCGA_GEO_results.Rdata         saved intermediates for visualize_results.R
2-GSEA/GO_ORA_*.csv, KEGG_ORA_*.csv  ORA result tables
2-GSEA/GO_GSEA_*.csv, KEGG_GSEA_*    GSEA result tables
3-Visualization/                     Volcano_Tumor_vs_Normal.pdf,
                                     Expression_boxplot_<gene>.pdf,
                                     GO/KEGG ORA + GSEA suites,
                                     ThemeEnrichment/ (theme dot-heatmaps, single-term gseaplot2)
6-Survival/                          clinical_clean.csv,
                                     KM_<gene>.pdf, KM_quartile_<gene>.pdf,
                                     KM_clinical_<var>.pdf, clinical_KM_summary.csv,
                                     univariate_Cox.csv/.pdf, multivariate_Cox.csv/.pdf
Analysis_summary.txt                 text summary of the run
sessionInfo.txt                      package versions (reproducibility)
```

### Notebook → numbered-layout mapping

The source notebook wrote everything to a flat `OUTDIR`. This runner routes:

- **DEG tables** (`DESeq2_Tumor_vs_Normal.csv`) → `1-DEG/`
- **ORA/GSEA tables** → `2-GSEA/`; their **figures** (suites, theme maps,
  single-term gseaplot2) → `3-Visualization/`
- **Expression boxplots + volcano** → `3-Visualization/`
- **All survival outputs** — both tables (`univariate_Cox.csv`,
  `multivariate_Cox.csv`, `clinical_KM_summary.csv`, `clinical_clean.csv`) and
  figures (`KM_*.pdf`, Cox forest PDFs) → `6-Survival/`. Keeping survival
  figures beside their tables also keeps `run_clinical_km()`'s co-written
  `KM_clinical_<var>.pdf` + `clinical_KM_summary.csv` together.

### Median-split KM de-duplication

The notebook computed the per-gene **median-split KM twice** — once in §7 and
again in §7.5. This runner keeps **one** median-split KM path: §7 loops
`GENES_FOR_SURVIVAL` a single time, producing the expression boxplot, the
median-split KM, and the quartile-split KM (`stratify_by_quantile`) together,
then §7.5's clinical KM, univariate Cox, and multivariate Cox follow. Genes
absent from the expression matrix are skipped with `warning(...); next`
(the notebook's tolerant behavior is preserved).

## Notes

- `SPECIES` selects the OrgDb (`org.Hs.eg.db` / `org.Mm.eg.db`) and KEGG organism
  code automatically — the notebook hard-coded human; here it is configurable.
- Expression for plotting/survival is `log2(TPM + 1)` when TPM is available,
  otherwise VST from the validated raw counts — never `log2(raw counts + 1)`,
  which leaves library size as a confounder.
- Survival runs only when the clinical table carries `vital_status`,
  `days_to_death`, and `days_to_last_follow_up`; otherwise that section is
  skipped with a message (e.g. a bare GEO phenotype table).
- Multivariate Cox needs ≥ 2 survival genes present in the expression matrix.
- The shared HTML report (`reports/analysis_report.qmd`) is oriented to the
  General layout, so `GENERATE_HTML_REPORT` defaults to `FALSE` here; the wiring
  is present if you add a TCGA/GEO-specific report template.

## Verifying against the demo

A runnable demo ships in `examples/demo_RNAseq_TCGA_GEO/`. Regenerate its local
inputs, point a config at them in local-file mode, and run:

```bash
Rscript examples/demo_RNAseq_TCGA_GEO/regenerate_demo_data.R
# then write a config with DOWNLOAD_FROM_GDC/GEO = FALSE and the demo paths,
# and: RNASEQ_LIB_DIR=<repo>/RNAseq_lib Rscript run_analysis.R /path/to/cfg.R
```
