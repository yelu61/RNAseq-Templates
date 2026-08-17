# RNAseq_lib Function Catalog

This catalog lists the helper functions in `RNAseq_lib/` grouped by file. Each entry gives the function signature and a one-line purpose. For detailed behavior, see the inline comments in each source file.

## `batch_utils.R` — batch effect diagnostics

- `plot_pca_by_batch_pdf(vsd, batch_vec, filename, ...)`
  PCA colored by batch vector; saves PDF.
- `summarize_pve_by_batch(vsd, batch_vec, condition_vec, ntop, npcs)`
  Returns variance-explained by batch (and condition) for top PCs.
- `plot_batch_pve_pdf(pve_df, filename, ...)`
  Barplot of batch/condition PVE across PCs.

## `design_utils.R` — paired/repeated-measures design

- `make_paired_col_data(sample_names, groups, group_levels, pair_id)`
  Build `colData` with `condition` and `pair_id`.
- `build_paired_design_formula(condition_col, pair_col)`
  Return `~ pair_id + condition`.
- `validate_paired_design(sample_names, groups, pair_id, group_levels)`
  Check pair completeness across conditions.

## `deg_utils.R` — DESeq2 DEG logic

- `validate_threshold_grid(threshold_grid, default_threshold)`
  Standardize and validate `THRESHOLD_GRID`.
- `run_deseq2_model(count_data, col_data, design_formula)`
  Build DESeqDataSet and run `DESeq()`.
- `extract_deseq2_results(dds, comparisons, ...)`
  Extract + shrink DESeq2 contrasts; returns list of result tables.
- `validate_pvalue_column(res_df, pvalue_column)`
  Check p-value column exists.
- `validate_lfc_column(res_df, lfc_column)`
  Check LFC column exists.
- `validate_numeric_column(res_df, column, label)`
  Check numeric column exists.
- `mark_deg_by_threshold(res_df, pvalue_thresh, log2fc_thresh, ...)`
  Add `significance` column (Up/Down/Not_Sig).
- `build_deg_threshold_sets(res_list, threshold_grid, ...)`
  Apply multiple thresholds.
- `summarize_deg_thresholds(deg_by_threshold, ...)`
  Count DEGs per threshold/comparison.
- `get_threshold_result(deg_by_threshold, threshold_name)`
  Extract results for one threshold.
- `write_all_gene_deg_results(res_list, outdir)`
  Write all-gene CSVs per comparison.
- `write_deg_threshold_outputs(deg_by_threshold, outdir, ...)`
  Write threshold-specific CSVs and summary.
- `summarize_lfc_strategies(res_list, threshold_grid, ...)`
  Compare raw vs shrunken LFC DEG counts.
- `write_lfc_strategy_summary(...)`
  Export raw-vs-shrunken summary.
- `summarize_deg_diagnostics(...)` / `write_deg_diagnostic_summary(...)`
  Comprehensive p-value/LFC diagnostic tables.
- `genes_for_enrichment(res_df, ...)`
  Return `sig`/`up`/`down` gene symbols.
- `ranked_gene_list(res_df, rank_column)`
  Named sorted vector for GSEA.
- `build_deg_gene_sets(...)`
  Build named gene sets for overlap.
- `prepare_upset_df(gene_sets)`
  Convert gene sets to UpSetR input.
- `plot_deg_upset_pdf(gene_sets, filename, ...)`
  UpSetR plot of DEG sets.
- `jaccard_overlap_matrix(gene_sets)`
  Pairwise Jaccard matrix.
- `plot_deg_overlap_heatmap_pdf(gene_sets, filename, ...)`
  Heatmap of Jaccard overlaps.
- `extract_intersection_genes(gene_sets, set_names, mode)`
  Extract common/union genes.

## `enrichment_utils.R` — enrichment execution

- `map_symbols_to_entrez(symbols, org_db)`
  Symbol → Entrez ID mapping with deduplication.
- `make_entrez_ranked_list(gene_list, org_db)`
  Convert ranked gene list to Entrez-ranked list.
- `run_go_ora(symbols, org_db, universe, ...)`
  GO over-representation analysis.
- `run_kegg_ora(symbols, org_db, universe, organism, ...)`
  KEGG over-representation analysis.
- `run_go_gsea(entrez_list, org_db, ...)`
  GO GSEA.
- `run_kegg_gsea(entrez_list, organism, ...)`
  KEGG GSEA.
- `run_threshold_ora(res_list, threshold_grid, org_db, ...)`
  Multi-threshold ORA loop across comparisons.
- `enrich_result_to_df(enrich_result)`
  Safe conversion to data frame.
- `build_multi_comparison_enrich_df(result_map, ...)`
  Bind enrichment results across comparisons.

## `geo_utils.R` — GEO public data

- `download_geo_series_matrix(geo_accession, destdir, getGPL)`
  Download GEO SeriesMatrix via `GEOquery`.
- `parse_geo_series_matrix(gse)`
  Extract expression, phenotype, and feature data.
- `prepare_geo_counts(expr, feature, ...)`
  Collapse probes to gene symbols.

## `io_utils.R` — input/output, QC, preprocessing

- `read_count_table(input_file, input_format)`
  Read TSV/CSV/Excel count table.
- `summarize_raw_input(rawcount, gene_name_col, ...)`
  QC summary of raw input.
- `detect_count_columns(rawcount, gene_name_col, ...)`
  Identify count columns.
- `validate_sample_design(...)`
  Validate sample/group/comparison consistency.
- `calculate_sample_qc(count_data, col_data, group_col)`
  Compute library size, detected genes, zero fraction, correlation.
- `flag_sample_qc(sample_qc, ...)`
  Flag samples by thresholds.
- `apply_sample_exclusion(count_data, col_data, sample_exclude)`
  Drop specified samples.
- `build_count_matrix(rawcount, gene_name_col, ...)`
  Validate counts, deduplicate symbols.
- `preview_count_matrix(count_data, gene_name_col, n)`
  Print preview.
- `make_col_data(count_data, sample_names, groups, group_levels)`
  Build `colData`.
- `filter_low_count_genes(count_data, groups, min_count)`
  Low-count filter.
- `calculate_filter_retention(raw, filtered)`
  Library retention per sample.
- `write_preprocessing_summary(summary_file, rows)`
  Write preprocessing summary CSV.
- `write_deg_excel(deg_by_threshold, outdir, filename)`
  Multi-threshold DEG Excel export.
- `write_all_genes_excel(res_list, outdir, filename)`
  All-gene results Excel export.

## `data_utils.R` — unified data loading, validation, and gene conversion

- `read_expression_matrix(file, gene_column, sample_ids, file_format)`
  Load CSV/TSV/Excel expression matrix; auto-detect gene column, coerce numeric, deduplicate rows by mean expression.
- `read_metadata(file, sample_column, required_columns, group_column, group_levels, time_column, time_levels)`
  Load metadata; rename duplicated sample columns, validate required columns, factorize group/time columns.
- `validate_samples_match(expr_samples, meta_samples, context, strict_order)`
  Check sample overlap and ordering between expression and metadata.
- `detect_expression_scale(mat, sample_n, gene_n)`
  Heuristic classification of expression matrix scale: raw_counts / log2_tpm / tpm / vst / unknown.
- `validate_expression_contract(mat, expected, tolerance)`
  Enforce declared expression units and gene/sample-name integrity; heuristic detection never chooses a transformation automatically.
- `counts_to_tpm(counts_mat, gene_lengths_kb)`
  Convert raw counts to TPM using gene lengths (kb).
- `counts_to_fpkm(counts_mat, gene_lengths_kb)`
  Convert raw counts to FPKM.
- `extract_gene_lengths(raw_annot, id_col, length_col, start_col, end_col, length_unit)`
  Extract gene lengths in kb from annotation columns.
- `collapse_by_symbol(symbols, mat, extra)`
  Collapse a count/expression matrix to unique symbols, keeping the highest-total-signal row per symbol; optionally carries a per-row annotation (e.g. gene length) through to the retained row for TPM/TME-ready tables.
- `validate_count_matrix(mat, require_integer, require_non_negative, min_samples)`
  Validate numeric, non-negative, optionally integer count matrix.
- `detect_gene_id_type(ids)`
  Heuristic: Ensembl IDs vs symbols.
- `convert_gene_ids(ids, from, to, species, method)`
  Convert gene IDs via AnnotationDbi / babelgene; supports human/mouse Ensembl↔symbol and MGI→HGNC.
- `convert_expression_rownames(expr, species, target, method)`
  Convert expression matrix row names and deduplicate by mean expression.

## `limma_voom_utils.R` — limma-voom DE

- `prepare_dge_for_voom(counts, group, ...)`
  Build DGEList, filter, normalize.
- `run_voom(dge, design, plot_file)`
  Run `limma::voom`.
- `remove_batch_effect_voom(v, batch, ...)`
  Visualization-only batch removal. Differential testing should include batch in `make_group_design()`.
  Apply `limma::removeBatchEffect`.
- `run_limma_contrasts(v, design, comparisons)`
  Fit + contrast + eBayes; standardized result list.
- `make_group_design(group)`
  No-intercept design matrix.
- `summarize_limma_deg(...)` / `write_limma_results(...)`
  Counts and CSV output.
- `ranked_gene_list_limma(res_df, rank_column)`
  Ranked list for GSEA.

## `plot_utils.R` — plotting & visualization

- `RNAseq_FIGURE_SPEC`
  Final-size journal dimensions, typography, raster resolution and font stack.
- `publication_dimensions(column, height_mm)` / `mm_to_in(x)`
  Convert final physical dimensions to export dimensions.
- `resolve_publication_font(preferred)`
  Select an installed Helvetica-compatible sans face with a portable fallback.
- `theme_publication(base_size, base_family)`
  Final-size 8 pt publication theme with a consistent typographic hierarchy.
- `make_group_colors(group_levels)`
  Generate group color palette.
- `save_pdf_plot(plot, filename, width, height)`
  Save ggplot atomically via cairo/pdf; reject NULL and near-empty device output.
- `save_pdf_device(filename, width, height, draw)`
  Transactional PDF export for grid/base graphics such as ComplexHeatmap,
  pheatmap, Mfuzz, UpSet, voom and survminer.
- `save_plot_bundle(plot, filename_stem, width_mm, height_mm, formats, dpi)`
  Export selected figures as PDF/SVG plus optional TIFF/PNG at final size.
- `wrap_term_labels(x, width, max_lines)`
  Wrap long labels at word boundaries, cap line count and add an ellipsis.
- `parse_ratio_numeric(x)`
  Parse "a/b" ratios.
- `zscore_rows(mat, cap)`
  Row Z-score with optional capping.
- `plot_pca_pdf(vsd, group_levels, group_colors, filename, intgroup)`
  PCA PDF.
- `plot_sample_distance_pdf(vsd, filename)`
  Sample distance heatmap.
- `plot_sample_qc_pdf(sample_qc, filename, group_colors)`
  Sample QC barplots.
- `plot_sample_correlation_pdf(sample_qc, filename)`
  Spearman correlation heatmap.
- `plot_count_distribution_pdf(count_data, filename)`
  Count distribution boxplot.
- `plot_filter_retention_pdf(retention_df, filename)`
  Filter retention barplot.
- `plot_deg_summary_pdf(deg_summary, filename)`
  DEG count summary.
- `plot_lfc_strategy_summary_pdf(lfc_strategy_summary, filename)`
  Raw vs shrunken LFC summary.
- `plot_volcano_pdf(res_df, comp_name, ...)`
  EnhancedVolcano PDF.
- `plot_expression_heatmap_pdf(mat, filename, ...)`
  ComplexHeatmap expression heatmap.
- `format_p_for_label(p)`
  Format p-value labels.
- `pairwise_effect_table(...)`
  Pairwise test summary.
- `plot_group_boxplot_pdf(..., stat_table = NULL)`
  Faceted boxplot with computed pairwise stats or an exact precomputed statistic table.
- `plot_group_violin_boxplot_pdf(..., stat_table = NULL)`
  Violin + boxplot with computed pairwise stats or an exact precomputed statistic table.
- `plot_group_bar_sem_pdf(...)`
  Mean ± SEM barplot with adjusted pairwise stats.
- `prepare_enrich_df(enrich_result, show_category)`
  Convert enrichment result to plot-ready df.
- `plot_enrich_dotplot(...)`
  Table-compatible enrichment dotplot with bounded external term labels.
- `plot_enrich_barplot_pdf(...)`
  Enrichment barplot.
- `plot_enrich_bidirectional_barplot_pdf(...)`
  Up/Down bidirectional barplot.
- `plot_enrich_network_pdf(...)`
  cnetplot + emapplot.
- `plot_enrich_suite_pdf(...)`
  Run dotplot + barplot + network suite.
- `plot_comparecluster_dotplot_pdf(...)`
  compareCluster dotplot.
- `plot_gsea_nes_barplot_pdf(...)`
  GSEA NES barplot.
- `plot_gsea_suite_pdf(...)`
  GSEA overview suite (dotplot, NES barplot, ridgeplot). Running curves are
  intentionally generated one term per file by `plot_gsea_term_figures_from_df()`.
- `plot_gsea_term_figure_pdf(...)`
  Single-term GSEA running-enrichment figure.
- `plot_gsea_term_figures_from_df(...)`
  Batch single-term figures.
- `default_enrichment_themes()` / `match_enrichment_themes()` / `prepare_theme_dotplot_df()`
  Theme dictionary and helpers.
- `plot_theme_dotheatmap_pdf(...)` / `plot_theme_dotheatmap_from_results(...)`
  Theme dot-heatmap with bounded automatic PDF dimensions.

## `pathway_utils.R` — gene-set scoring & pathway visualization

- `score_gene_sets(expr, gene_sets, method, kcdf, min_size, verbose)`
  Validate and score custom gene sets via GSVA; defaults to five overlapping features and attaches an overlap audit.
- `pathway_group_comparison(score_mat, group, group_levels, comparisons, ...)`
  Sample-name-aligned per-set group comparison; delta = treatment - control, global BH.
- `read_gene_set_registry(path, ...)` / `audit_gene_set_registry(gene_sets, ...)`
  Preserve version/source metadata and audit registered-gene overlap with the analyzed expression matrix.
- `plot_pathway_delta_summary_pdf(comp_df, filename, ...)`
  Bidirectional delta bar chart across pathways x comparisons with significance stars.
- `plot_pathway_sensitivity_matrix_pdf(comp_df, filename, ...)`
  Complete module-by-contrast matrix combining delta color with parametric and exact-permutation adjusted-P labels.
- `plot_gene_set_registry_qc_pdf(audit_df, filename, ...)`
  Registry coverage chart showing expressed-gene overlap fraction, overlap/input denominators, and source type.
- `plot_keygenes_log2fc_heatmap_pdf(fc_long, filename, row_group, ...)`
  Gene x comparison log2FC heatmap with preserved NA values, duplicate-key validation, and optional row grouping.
- `melt_gene_expression(expr, genes, group)`
  Long-format selected genes for grouped plots.
- `plot_gene_expression_pdf(expr, genes, group, group_levels, group_colors, filename, ...)`
  Per-gene grouped expression plots (faceted boxplot or per-gene box/violin).

## `tme_utils.R` — TME deconvolution helpers

- `undo_log_expr(expr, is_log, log_base)`
  Reverse log transformation for TME tools that expect non-log input.
- `validate_tme_input(expr)`
  Validate expression matrix: numeric, non-negative, rownames present; warn on Ensembl IDs.
- `deduplicate_expression_by_symbol(expr, symbol_col)`
  Deduplicate rows by gene symbol, keeping the row with highest mean expression.
- `convert_mouse_symbols_to_human(expr, mouse_attr, human_attr, host, fallback_hosts, ortholog_cache, verbose)`
  Convert mouse MGI symbols to human HGNC symbols via `biomaRt::getLDS`. Tries
  `host` then `fallback_hosts`; with `ortholog_cache` (`.rds` path) the mapping is
  cached so repeat runs are offline-deterministic (only uncached genes hit network).
- `prepare_tme_expression(expr, is_log, species, log_base, ortholog_cache, host, fallback_hosts, verbose)`
  One-stop preparation: undo log + optional mouse-to-human conversion. Set
  `ortholog_cache` for offline-deterministic repeat runs of mouse data.
- `calc_tme_barplot_size(n_samples, ...)` / `calc_tme_boxplot_size(n_celltypes, ...)`
  Automatic figure-size calculation based on sample/cell-type counts.
- `read_cibersort_signature(signature_file)`
  Read a CIBERSORT signature file (LM22.txt or cibersort_mouse_22.csv) into a numeric matrix.
- `run_native_cibersort(expr, signature_file, cibersort_script, is_log, perm, QN, id_column, verbose)`
  Run the bundled standalone CIBERSORT implementation on a prepared expression matrix. Supports both human LM22 and the mouse `cibersort_mouse_22.csv` signature.
- `run_iobr_deconvolution(expr, methods, perm, arrays, id_column)`
  Run IOBR methods (`cibersort`, `epic`, `xcell`, `estimate`).
- `combine_tme_results(result_list, id_column)`
  Inner-join IOBR result tables by sample ID.
- `melt_tme_results(tme_df, id_column, group_df, sample_col, group_col)`
  Convert wide IOBR result to long format for plotting.
- `melt_estimate_scores(estimate_df, id_column, group_df, sample_col, group_col)`
  Convert ESTIMATE scores to long format.
- `compare_native_iobr_cibersort(native_df, iobr_df, id_column, method)`
  Align native and IOBR CIBERSORT results, compute per-cell-type correlation/RMSE/MAE/paired t-test, and return summary + long-format data.
- `plot_cibersort_correlation_pdf(long_df, filename, title, width, height)`
  Faceted scatter plot of native vs IOBR fractions with per-cell-type Pearson correlation.
- `plot_cibersort_difference_pdf(long_df, filename, title, width, height)`
  Bland–Altman-style difference plot (native − IOBR vs mean) per cell type.
- `plot_tme_barplot_pdf(long_df, group_col, sample_col, filename, ...)`
  Stacked barplot of TME cell fractions/scores.
- `plot_tme_boxplot_pdf(long_df, group_col, value_col, filename, ...)`
  Faceted box/violin plots per cell type with pairwise statistics.
- `plot_tme_per_celltype_pdf(long_df, group_col, value_col, filename_prefix, ...)`
  One PDF per cell type for focused single-cell-type figures.
- `plot_tme_heatmap_pdf(tme_df, meta, group_col, sample_col, filename, ...)`
  pheatmap wrapper for xCell / ESTIMATE scores with group annotation colors.
- `plot_estimate_boxplot_pdf(long_df, group_col, filename, ...)`
  ESTIMATE score boxplots by group.

## `survival_utils.R` — survival analysis

- `validate_surv_df(surv_df, time_col, status_col)`
  Validate survival data frame.
- `run_univariate_cox(surv_df, vars, ...)`
  Per-variable Cox models.
- `run_multivariate_cox(surv_df, vars, ...)`
  Multivariate Cox model.
- `plot_cox_forest_pdf(cox_df, filename, ...)`
  Cox forest plot.
- `plot_km_by_group_pdf(surv_df, group_col, filename, ...)`
  Kaplan–Meier by categorical group.
- `stratify_by_quantile(x, n_groups, labels)`
  Quantile-based grouping.

## `tcga_utils.R` — TCGA public data

- `build_tcga_query(project, ...)`
  Build `TCGAbiolinks::GDCquery`.
- `extract_tcga_assays(se, assay_names)`
  Pull assay matrices.
- `symbolize_and_dedup(mat, id_map, ...)`
  Ensembl → symbol, dedup by mean expression.
- `build_id_map_from_se(se, ...)`
  Derive id→symbol map from `rowData(se)`.
- `extract_tcga_clinical(se)`
  Extract and flatten clinical `colData`.
- `infer_tcga_tumor_normal(barcodes, ...)`
  Infer tumor/normal from TCGA barcode.
- `prepare_tcga_survival(clinical, ...)`
  Build survival data frame.
- `run_clinical_km(surv_df, clinical_df, var_cols, ...)`
  KM for clinical variables.
- `plot_km_by_median_pdf(surv_df, value_col, filename, ...)`
  KM by median split.
- `plot_tcga_gene_boxplot_pdf(expr_df, gene, group_vec, filename, ...)`
  Single-gene expression boxplot.

## `timecourse_utils.R` — time-course / temporal analysis

- `aggregate_expr_by_group(expr, group_vec, fun)`
  Aggregate expression by group means.
- `prepare_mfuzz_eset(expr, ...)`
  Prepare Mfuzz ExpressionSet.
- `run_mfuzz(eset, n_clusters, seed, m)`
  Run Mfuzz soft clustering.
- `extract_mfuzz_clusters(mfuzz_result, eset, min_acore)`
  Cluster assignments and membership.
- `plot_mfuzz_trends_pdf(eset, mfuzz_result, filename, ...)`
  Cluster trend PDF.
- `plot_timecourse_heatmap_pdf(expr, cluster_df, ...)`
  Z-score heatmap split by cluster.
- `run_mfuzz_cluster_ora(cluster_df, org_db, ...)`
  ORA per cluster.
- `run_timepoint_vs_baseline_deseq2(count_data, col_data, ...)`
  DESeq2 time-point vs baseline contrasts.
- `summarize_timepoint_deg(res_list, ...)`
  DEG counts per time point.
- `write_timepoint_deg_results(res_list, outdir, ...)`
  Export time-point DEG tables.
- `plot_timepoint_deg_summary_pdf(summary_df, filename)`
  Barplot of time-point DEG counts.
- `write_mfuzz_cluster_table(cluster_df, filename)`
  Write cluster assignments.
- `summarize_mfuzz_clusters(cluster_df)`
  Cluster size summary.

## `report_utils.R` — unified HTML report

- `render_analysis_report(outdir, report_file, template, params)`
  Render `reports/analysis_report.qmd` into a single self-contained HTML report from the saved DEG/ORA/GSEA/QC outputs (no re-analysis). A `.qmd` uses Quarto when available; otherwise it automatically falls back to the sibling `reports/analysis_report.Rmd` twin via `rmarkdown` (same body, rmarkdown YAML).
- `list_report_figures(dir, pattern, base)`
  List PDF figures under a directory as portable relative paths.
- `read_csv_safe(path, ...)`
  Read a CSV if it exists, else return `NULL` (for conditional report tables).
- `default_report_coverage_rules()` / `build_report_coverage_manifest(...)`
  Define and audit the required/optional analysis domains represented in a comprehensive report.
- `write_report_coverage_manifest(...)` / `validate_report_coverage(...)`
  Persist the coverage contract and reject unresolved required domains or undocumented omissions.
- `validate_claim_evidence_ledger(...)`
  Validate project-specific conclusions, evidence levels, alternatives, falsifiers, next evidence, and source-file provenance before report rendering.

## TME constants

- `immune_gene_sets`
  Built-in 28-cell-type immune signature collection (Charoentong et al. 2017, human symbols) used as the default ssGSEA gene sets in the TME deconvolution template.
