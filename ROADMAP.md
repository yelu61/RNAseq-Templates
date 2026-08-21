# RNAseq-Templates Roadmap

## Status

This project provides production-oriented bulk RNA-seq analysis templates with
a shared R helper library (`RNAseq_lib/`). The command-line runner is preferred
for reproducible project execution; notebooks remain available for interactive
exploration. Every template — General plus the five topic templates — now ships
a `config.R` + `run_analysis.R` + `visualize_results.R` runner under
`templates/`, the core pipeline and six topic templates are functional, and all
helper modules have unit tests. `tools/notebook_to_runner.R` converts a
notebook's parameter cell into a `config.R` + runner draft.

**Current release: v0.9.0** (see [CHANGELOG.md](CHANGELOG.md)).

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
- [x] Added final-size publication typography, safe multi-format export, and
      visual QA rules across shared plotting helpers
- [x] Added a two-layer real-project output contract based on the
      `202607XXR_RPE` implementation audit
- [x] Added separate complete-execution and read-only run-review General notebooks
- [x] Added run manifests, duplicate registry and explicit retention lifecycle
- [x] Added independent report re-rendering plus scientific-review publish gate

## Planned / Future

- [x] Add automated CI/GitHub Actions to run `examples/run_demo_smoke_test.R` on push
- [x] Add formal unit tests for all `RNAseq_lib` helpers
- [x] Command-line runner for the General template (`templates/General/run_analysis.R` + `config.R`)
- [x] Optional TME deconvolution switch (`RUN_TME`) inside `run_analysis.R` (output to `4-TME/`)
- [x] Targeted re-visualization script (`templates/General/visualize_results.R`) from saved results
- [x] Fix paired-design and batch-covariate integration bugs in the General template
- [x] Get GitHub Actions smoke test fully green (dependency install + all six template demos)
- [x] Add limma-voom and WGCNA demos to the smoke test (now exercises all six templates)
- [x] Create a converter to turn notebook parameter cells into `config.R` + a
      `run_analysis.R` draft for headless execution (`tools/notebook_to_runner.R`,
      R + jsonlite; emits a conversion report for manual refinement)
- [ ] Add more public-data examples (GEO local-file mode, TCGA local-file mode)
- [x] Add interactive HTML report generation (`reports/analysis_report.qmd`, rendered via quarto or rmarkdown)
- [ ] Add single-cell reference integration / deconvolution (e.g., Bisque, CIBERSORTx)
- [ ] Support additional model organisms (rat, zebrafish) via configurable `org.*.eg.db`
- [ ] Introduce `renv.lock` for long-term dependency reproducibility (deferred; platform-specific binaries and GitHub/R-Forge sources make this non-trivial)
- [ ] Publish v1.0.0
- [x] Provide example datasets for each topic template (General, limma-voom, WGCNA, TME, TimeCourse, TCGA/GEO all done)
- [x] Add runnable demos for RNAseq_TME_Deconvolution_Template.ipynb
- [x] Add runnable demos for RNAseq_TimeCourse_Template.ipynb
- [x] Add runnable demos for RNAseq_TCGA_GEO_Template.ipynb (local-file mode)
- [x] Add unified HTML report generation (`reports/analysis_report.qmd` + `render_analysis_report()`)
- [x] Add command-line runners (`config.R` + `run_analysis.R` +
      `visualize_results.R`) for all five topic templates (Limma_Voom, WGCNA,
      TME, TimeCourse, TCGA_GEO), matching the General runner conventions;
      topic demos now drive these runners so CI exercises them

## How to contribute

1. Open an issue or discussion describing the improvement.
2. Prefer adding helpers to `RNAseq_lib/` rather than copying code into notebooks.
3. Update `references/FUNCTION_CATALOG.md` when adding new helpers.
4. Add or update `references/TROUBLESHOOTING.md` for common user issues.
5. If a new dependency is required, add it to `install_dependencies.R`.
