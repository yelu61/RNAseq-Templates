# Changelog

All notable changes to RNAseq-Templates are documented in this file.

## [Unreleased]

### Added

- `examples/run_demo_smoke_test.R`: automated smoke test that runs the General notebook core pipeline on demo data and validates outputs.
- `RNAseq_lib/batch_utils.R`: batch-effect PCA (`plot_pca_by_batch_pdf`), batch variance-explained summary (`summarize_pve_by_batch`), and PVE barplot (`plot_batch_pve_pdf`).
- `RNAseq_lib/design_utils.R`: paired/repeated-measures design helpers (`make_paired_col_data`, `build_paired_design_formula`, `validate_paired_design`).
- `RNAseq_lib/geo_utils.R`: GEO SeriesMatrix download and parsing helpers (`download_geo_series_matrix`, `parse_geo_series_matrix`, `prepare_geo_counts`).
- `RNAseq_lib/io_utils.R`: Excel export helpers (`write_deg_excel`, `write_all_genes_excel`).
- `RNAseq_lib/timecourse_utils.R`: time-point vs baseline DESeq2 helpers (`run_timepoint_vs_baseline_deseq2`, `summarize_timepoint_deg`, `write_timepoint_deg_results`, `plot_timepoint_deg_summary_pdf`).
- `RNAseq_lib/tcga_utils.R`: `build_id_map_from_se()` to derive ENSEMBL→symbol maps from `SummarizedExperiment::rowData`, and `run_clinical_km()` for clinical-variable Kaplan–Meier curves.
- `RNAseq_PALETTE` constant and visualization style guide at `references/VISUALIZATION_STYLE_GUIDE.md`.
- Reference documentation: `references/PARAMETER_REFERENCE.md`, `references/FUNCTION_CATALOG.md`, `references/TEMPLATE_SELECTION.md`, `references/TROUBLESHOOTING.md`.
- `RNAseq_lib/timecourse_utils.R`: restored orphaned code fragment as `write_mfuzz_cluster_table()`.
- `examples/demo_data/demo_counts.tsv`: regenerated with real mouse gene symbols so enrichment steps pass in the smoke test.
- `examples/run_demo_smoke_test.R`: made enrichment assertions data-dependent; now verifies that ORA summary exists rather than requiring a specific threshold-specific CSV.
- `.github/workflows/smoke-test.yml`: GitHub Actions workflow to run the smoke test and unit tests on push/PR.
- `tests/testthat/` and `tests/testthat.R`: initial unit-test skeleton covering deg, io, and enrichment helpers.
- `examples/regenerate_demo_data.R`: helper to regenerate demo count table with real mouse symbols.
- `examples/demo_RNAseq_limma_voom/` and `examples/demo_RNAseq_WGCNA/`: pre-configured runnable demos for the limma-voom and WGCNA templates.

### Changed

- `RNAseq_lib/plot_utils.R`: `plot_deg_summary_pdf()` now tolerates single-threshold summaries (no `Threshold` column).
- `notebooks/RNAseq_WGCNA_Template.ipynb`: force first column to character when `GENE_COLUMN` is NULL, and convert expression matrix to numeric matrix before WGCNA input, preventing failures with numeric-looking gene symbols.

- `RNAseq_lib/deg_utils.R`: fixed `extract_deseq2_results()` referencing `raw_df` before it was defined.
- `notebooks/RNAseq_limma_voom_Template.ipynb`: replaced broken `BATCH_COLUMN` logic with explicit `BATCH_VECTOR` parameter.
- `notebooks/RNAseq_General.ipynb`:
  - Added `BATCH_VECTOR`, `PAIR_ID`, and `EXPORT_EXCEL` parameters.
  - Sources `batch_utils.R` and `design_utils.R`.
  - Supports paired design when `PAIR_ID` is supplied.
  - Generates batch-effect diagnostics when `BATCH_VECTOR` is supplied.
  - Exports `1-DEG/DEG_results.xlsx` when `EXPORT_EXCEL = TRUE`.
- `notebooks/RNAseq_TimeCourse_Template.ipynb`:
  - Added raw-count and metadata parameters for Section 6.
  - Implemented runnable time-point vs baseline DEG with optional paired design.
- `notebooks/RNAseq_TCGA_GEO_Template.ipynb`:
  - Added `DOWNLOAD_FROM_GEO` / `GEO_ACCESSION` parameters.
  - `GENE_ID_MAP_FILE` now defaults to `NULL` and falls back to `build_id_map_from_se()`.
  - Added `CLINICAL_VARS_FOR_KM` parameter and clinical-variable KM output.
- `install_dependencies.R`: added `ashr`, `ggrepel`, `data.table`, `matrixStats`, `msigdbr`, `DOSE`, `BiocParallel`, and `GEOquery`.
- `RNAseq_lib/enrichment_utils.R` and `RNAseq_lib/plot_utils.R`: added friendly `message()` logs when enrichment or plots are skipped due to empty results.

### Fixed

- `extract_deseq2_results()` no longer fails with "object 'raw_df' not found".
- limma-voom batch correction no longer checks nonsensical `names(GROUPS)` / `Sys.getenv()` conditions.

## Earlier releases

- Initial public release with six notebook templates, `RNAseq_lib` helpers, demo data, and Chinese Getting Started guide.
- Added IOBR TME, TCGA/GEO, limma-voom, time-course, and survival analysis modules.
- Added publication-grade enrichment visualization (theme dot-heatmaps and single-term GSEA figures).
