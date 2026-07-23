# Changelog

All notable changes to RNAseq-Templates are documented in this file.

## [Unreleased]

### Added

- **Command-line runner for the General template** under `templates/General/`:
  - `run_analysis.R` + `config.R` run the full General pipeline (all 15 sections) non-interactively via `Rscript`, driven by a single editable config file. `RNAseq_lib` is located via `RNASEQ_LIB_DIR`, the project directory, the repository root, or the parent directory. Intended for reproducible / batch / headless execution alongside the notebook.
  - **Optional TME deconvolution switch** (`RUN_TME`) in `run_analysis.R`: builds TPM from raw counts + gene lengths, then runs native ESTIMATE, IOBR (`estimate`/`cibersort`/`epic`/`xcell`), and ssGSEA immune signatures, writing to `4-TME/`. Sub-switches `RUN_TME_ESTIMATE` / `RUN_TME_IOBR` / `RUN_TME_SSGSEA`; IOBR/estimate absence degrades gracefully. Mouse input is converted to human orthologs via biomaRt.
  - `run_analysis.R` now caches `gseaResult` objects to `2-GSEA/gsea_results.rds` and saves `GROUP_LEVELS`/`COMPARISONS`/`SPECIES`/`colData` into `1-DEG/DEG_results.Rdata` so downstream figures can be regenerated without recompute.
  - **`visualize_results.R`**: targeted re-visualization from saved results with no DESeq2/ORA/GSEA recompute. Three independent sections — key-gene bar/SEM + heatmap (from `DEG_results.Rdata`), single-term gseaplot2 figures (from cached `gsea_results.rds`), and ORA theme dot-heatmaps (from saved ORA csvs).
  - `templates/General/README.md` documenting usage, batch execution, switches, and when to use the notebook instead.
- **Unit tests for `plot_utils.R`** (`tests/testthat/test-plot_utils.R`, 23 test blocks): palette/theme, label/ratio/z-score/SEM helpers, `pairwise_effect_table`, `prepare_enrich_df`, and theme matching.

### Fixed

- `plot_gsea_nes_barplot_pdf()` no longer fails with `factor level is duplicated` when KEGG GSEA results carry `NA` or duplicated `Description` values (the online KEGG map supplies IDs without names); missing labels now fall back to the term ID and are deduplicated.
- `build_multi_comparison_enrich_df()` no longer drops enrichment rows whose `ONTOLOGY` is `NA` when an `ontology_filter` is set — this previously blanked ORA theme dot-heatmaps built from csv-read results (e.g. in `visualize_results.R`).
- Native-ESTIMATE table parsing in `run_analysis.R` now selects real sample columns explicitly instead of positionally, so the extra `Description.1` column ESTIMATE writes is not mistaken for a sample.
- Removed a duplicated `tme_utils.R` section in `references/FUNCTION_CATALOG.md`.
- `Rplots.pdf` added to `.gitignore` (stray R plotting side-effect).

### Added

- **Runnable demos for all remaining topic templates**, each with a `regenerate_demo_data.R` + `run_demo.R` pair under `examples/`:
  - `examples/demo_RNAseq_TME/`: raw counts + gene lengths → TPM → ESTIMATE / IOBR (estimate, cibersort, epic) / ssGSEA. Uses human symbols so the whole run is offline.
  - `examples/demo_RNAseq_TimeCourse/`: VST expression across 4 time points with injected temporal patterns → Mfuzz clustering, plus raw-count time-point-vs-baseline DESeq2.
  - `examples/demo_RNAseq_TCGA_GEO/`: local-file mode (no network) with a synthetic TCGA-like cohort (TCGA barcodes, raw counts + TPM + clinical) → Tumor-vs-Normal DESeq2, KM survival, clinical KM, Cox, ORA/GSEA.
- `examples/run_demo_smoke_test.R` now drives every topic demo after the General demo and fails if any one fails, so each push exercises all templates.
- **Unified HTML report**: `reports/analysis_report.qmd` + `RNAseq_lib/report_utils.R` (`render_analysis_report()`). Assembles the saved DEG/ORA/GSEA/QC CSV and PDF outputs into a single self-contained HTML document without re-running the analysis. `RNAseq_General.ipynb` gains a `GENERATE_HTML_REPORT` parameter and a report cell (Section 15). Renders with the quarto CLI when available, else falls back to `rmarkdown`.
- `immune_gene_sets` built-in 28-cell-type immune signature collection (Charoentong et al. 2017) in `RNAseq_lib/tme_utils.R`, used by the TME ssGSEA step.
- CI (`smoke-test.yml`) installs the extra packages the new demos need (Mfuzz, WGCNA, survival, survminer, estimate, corrplot, e1071, babelgene, rprojroot) and rmarkdown for the report.

### Fixed

- **`immune_gene_sets` was undefined** in `RNAseq_TME_Deconvolution_Template.ipynb` — the ssGSEA cell referenced a variable that was never created and is not exported by IOBR, which would error for every user. Now provided as a built-in constant in `tme_utils.R`.
- **Native ESTIMATE failed with "找不到对象 'common_genes'"** in the TME template: `requireNamespace("estimate")` does not resolve the package's lazy-data objects, but `filterCommonGenes()`/`estimateScore()` reference them directly. Both the notebook and the demo now call `utils::data("common_genes"/"SI_geneset", package = "estimate")` first.
- **`plot_km_by_group_pdf()` failed with "object of type 'symbol' is not subsettable"** under R ≥ 4.x: the survfit call captured the formula as a local symbol, which `ggsurvplot()` could not re-evaluate. The formula is now inlined into the call via `eval(substitute(...))`. This broke all clinical-variable KM plots in the TCGA-GEO template.
- TimeCourse time-point-vs-baseline DEG no longer includes a single-level `condition` column in the DESeq2 design (which errored with "design contains variables with all samples having the same value"); the condition covariate is only added when it actually varies.

### Added

- **New `RNAseq_lib/data_utils.R`** with reusable data loading, validation, and gene-conversion helpers:
  - `read_expression_matrix()`, `read_metadata()`, `validate_samples_match()`
  - `detect_expression_scale()` for heuristic classification of raw counts / TPM / log-scale / VST
  - `counts_to_tpm()`, `counts_to_fpkm()`, `extract_gene_lengths()`, `validate_count_matrix()`
  - `convert_gene_ids()`, `convert_expression_rownames()` supporting human/mouse Ensembl↔symbol and MGI→HGNC
- **Unit tests** for `data_utils.R` in `tests/testthat/test-data_utils.R`.
- `validate_expression_contract()` enforces declared TPM/log2(TPM+1)/VST input units and gene/sample-name integrity.
- `tests/validate_notebooks.R` parses every notebook code cell and every shared R script; CI runs it on push and pull requests.
- `RNAseq_TME_Deconvolution_Template.ipynb` now defaults to raw integer counts plus gene lengths and computes TPM internally.
- `RNAseq_limma_voom_Template.ipynb`: new `SPECIES` parameter (default `"human"`) drives species-aware org.db and KEGG organism code selection.
- `RNAseq_TCGA_GEO_Template.ipynb`: new `GDC_COUNTS_ASSAY` and `GDC_TPM_ASSAY` parameters to explicitly select assays instead of relying on hard-coded `assays_list[[4]]`.

### Changed

- `RNAseq_TME_Deconvolution_Template.ipynb`, `RNAseq_TimeCourse_Template.ipynb`, `RNAseq_WGCNA_Template.ipynb`, and `RNAseq_TCGA_GEO_Template.ipynb` now load expression and metadata via `data_utils.R` helpers and perform explicit sample-name validation.
- `RNAseq_limma_voom_Template.ipynb` uses `validate_count_matrix()` after `build_count_matrix()`.
- limma-voom batch adjustment is now fitted as a model covariate; `removeBatchEffect()` is reserved for visualization.
- `RNAseq_General.ipynb` sources `data_utils.R` and calls `validate_count_matrix()` after building the count matrix.

### Fixed

- TME no longer treats VST/rlog as invertible log2(TPM+1); VST is rejected for CIBERSORT/EPIC/ESTIMATE input.
- Restored the truncated TME visualization cell and its CIBERSORT/EPIC/xCell/ESTIMATE outputs.
- Fixed successful gene-ID conversion/deduplication in `convert_expression_rownames()` and rejected unsafe numeric first-column inference.
- GEO SeriesMatrix assays must pass raw-integer-count validation before DESeq2; normalized GEO data are directed to limma.
- TPM conversion now rejects zero-total-RPK samples, and metadata rejects unknown factor levels.

- `plot_tme_heatmap_pdf()` now orders samples by group, then clusters within each group, so heatmaps keep biological replicates together while preserving within-group structure.
- `melt_estimate_scores()` now recognizes IOBR's `_estimate`-suffixed score columns and strips the suffix for consistent plotting.
- `get_cibersort_category_map("human")` now matches both canonical spaced LM22 names and IOBR's underscore-separated column names, fixing empty/wrong broad-category aggregation for IOBR CIBERSORT.
- `plot_estimate_boxplot_pdf()` supports `ncol`, `save_individual`, and `individual_prefix` for a combined 1×4 layout plus per-score single plots.

- **Mouse TME deconvolution support** in `RNAseq_TME_Deconvolution_Template.ipynb` and `RNAseq_lib/tme_utils.R`:
  - New `SPECIES` parameter (`"human"` / `"mouse"`) in the TME notebook.
  - `convert_mouse_symbols_to_human()` uses `biomaRt::getLDS` to map MGI symbols to HGNC symbols.
  - `prepare_tme_expression()` unifies log reversal + optional mouse-to-human conversion.
  - `deduplicate_expression_by_symbol()` keeps the highest-mean duplicate after ortholog conversion.
  - `validate_tme_input()` checks numeric/non-negative values, rownames, and warns on Ensembl IDs.
  - All TME methods (native ESTIMATE, native CIBERSORT, IOBR `estimate`/`cibersort`/`epic`/`xcell`, and ssGSEA immune signatures) now route through the prepared human-symbol matrix.
- **TME visualization improvements**:
  - New `GROUP_COLORS` parameter for user-defined group colors; falls back to `make_group_colors()`.
  - `calc_tme_barplot_size()` and `calc_tme_boxplot_size()` auto-adjust figure size by sample/cell-type counts.
  - New `plot_tme_per_celltype_pdf()` generates one focused violin/box PDF per cell type for CIBERSORT and EPIC.
  - New `plot_tme_heatmap_pdf()` wrapper produces xCell and IOBR ESTIMATE heatmaps with consistent group annotation colors.
  - All boxplots (ESTIMATE, CIBERSORT, EPIC, ssGSEA) now apply `group_colors` consistently.
- Metadata loader now renames duplicated `sample` columns in `colData.csv` to avoid downstream subsetting errors.
- `tests/testthat/test-tme_utils.R`: unit tests for TME helpers.
- `install_dependencies.R`: added `biomaRt` for mouse-to-human ortholog lookup.
- `references/PARAMETER_REFERENCE.md`: documented TME parameters and per-method input requirements.
- `references/FUNCTION_CATALOG.md`: documented new TME helpers.

### Changed

- `notebooks/RNAseq_TME_Deconvolution_Template.ipynb`:
  - Replaced per-method `undo_log_expr()` calls with a single `prepare_tme_expression()` step.
  - Native ESTIMATE, IOBR, native CIBERSORT, and ssGSEA now all use the prepared `expr_tme` / `expr_for_ssgsea` matrix.
  - Visualization block now uses dynamic sizing, per-cell-type outputs, and IOBR ESTIMATE heatmap.

### Fixed

- TME template no longer silently assumes human gene symbols; mouse data is now explicitly converted before deconvolution.
- TME plots no longer ignore user color preferences; `GROUP_COLORS` is propagated throughout.
- `plot_tme_heatmap_pdf()` now orders samples by group (using the factor levels) and clusters within each group, so heatmaps keep biological replicates together while preserving within-group structure.
- `melt_estimate_scores()` now recognizes IOBR's `_estimate`-suffixed score columns and strips the suffix for consistent plotting.
- `get_cibersort_category_map("human")` now matches both canonical spaced LM22 names and IOBR's underscore-separated column names, fixing empty/wrong broad-category aggregation for IOBR CIBERSORT.
- `plot_estimate_boxplot_pdf()` supports `ncol`, `save_individual`, and `individual_prefix` for a combined 1×4 layout plus per-score single plots.

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
