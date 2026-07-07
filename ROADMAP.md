# RNAseq-Templates Roadmap

## Status

This project provides notebook-first bulk RNA-seq analysis templates with a shared R helper library (`RNAseq_lib/`). The core pipeline and six topic templates are functional; ongoing work focuses on stability, documentation, and additional workflow helpers.

## Completed

- [x] General DESeq2 workflow with multi-threshold DEG (strict/standard/loose)
- [x] Publication-grade visualization theme and PDF outputs
- [x] ORA (GO/KEGG) and GSEA with theme dot-heatmaps and single-term figures
- [x] GSVA, single-gene plots, DEG overlap (UpSet + Jaccard), optional TF analysis
- [x] TME deconvolution (IOBR + ESTIMATE + ssGSEA)
- [x] WGCNA co-expression network template
- [x] TCGA/GEO public data mining template with survival KM/Cox
- [x] limma-voom alternative DE template
- [x] Time-course Mfuzz clustering template
- [x] Zero-programming-experience getting-started guide (Chinese)

## Recently implemented

- [x] Fixed `extract_deseq2_results()` undefined `raw_df` bug
- [x] Fixed limma-voom `BATCH_COLUMN` logic; replaced with `BATCH_VECTOR`
- [x] Synchronized `install_dependencies.R` with packages actually used in notebooks
- [x] Added demo smoke-test script: `examples/run_demo_smoke_test.R`
- [x] Improved empty-result logging in enrichment and plotting helpers
- [x] Implemented TimeCourse Section 6: time-point vs baseline DESeq2
- [x] Added batch-effect diagnostics (`batch_utils.R`)
- [x] Added paired-design helpers (`design_utils.R`)
- [x] Added multi-threshold DEG Excel export (`write_deg_excel`)
- [x] Added GEO SeriesMatrix download helpers (`geo_utils.R`)
- [x] Added TCGA id-map fallback from `rowData(se)` and clinical-variable KM
- [x] Created `references/` documentation and visualization style guide

## Planned / Future

- [x] Add automated CI/GitHub Actions to run `examples/run_demo_smoke_test.R` on push
- [x] Add formal unit tests for `RNAseq_lib` helpers
- [ ] Create a Python/CLI wrapper to convert notebook parameter cells to R scripts for headless execution
- [ ] Add more public-data examples (GEO local-file mode, TCGA local-file mode)
- [ ] Add interactive HTML report generation (using R Markdown or Quarto)
- [ ] Add single-cell reference integration / deconvolution (e.g., Bisque, CIBERSORTx)
- [ ] Support additional model organisms (rat, zebrafish) via configurable `org.*.eg.db`
- [x] Provide example datasets for each topic template (General, limma-voom, WGCNA done; TME, TimeCourse, TCGA/GEO planned)
- [ ] Add runnable demos for RNAseq_TME_Deconvolution_Template.ipynb
- [ ] Add runnable demos for RNAseq_TimeCourse_Template.ipynb
- [ ] Add runnable demos for RNAseq_TCGA_GEO_Template.ipynb (local-file mode)

## How to contribute

1. Open an issue or discussion describing the improvement.
2. Prefer adding helpers to `RNAseq_lib/` rather than copying code into notebooks.
3. Update `references/FUNCTION_CATALOG.md` when adding new helpers.
4. Add or update `references/TROUBLESHOOTING.md` for common user issues.
5. If a new dependency is required, add it to `install_dependencies.R`.
