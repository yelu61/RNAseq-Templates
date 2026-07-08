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

## `limma_voom_utils.R` — limma-voom DE

- `prepare_dge_for_voom(counts, group, ...)`  
  Build DGEList, filter, normalize.
- `run_voom(dge, design, plot_file)`  
  Run `limma::voom`.
- `remove_batch_effect_voom(v, batch, ...)`  
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

- `theme_publication(base_size, base_family)`  
  Default publication theme.
- `make_group_colors(group_levels)`  
  Generate group color palette.
- `save_pdf_plot(plot, filename, width, height)`  
  Save ggplot as PDF via cairo/pdf.
- `wrap_term_labels(x, width)`  
  Wrap long labels.
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
- `plot_group_boxplot_pdf(...)`  
  Faceted boxplot with pairwise stats.
- `plot_group_violin_boxplot_pdf(...)`  
  Violin + boxplot with pairwise stats.
- `plot_group_bar_sem_pdf(...)`  
  Mean ± SEM barplot with pairwise stats.
- `prepare_enrich_df(enrich_result, show_category)`  
  Convert enrichment result to plot-ready df.
- `plot_enrich_dotplot(...)`  
  Enrichment dotplot.
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
  GSEA plot suite.
- `plot_gsea_term_figure_pdf(...)`  
  Single-term GSEA running-enrichment figure.
- `plot_gsea_term_figures_from_df(...)`  
  Batch single-term figures.
- `default_enrichment_themes()` / `match_enrichment_themes()` / `prepare_theme_dotplot_df()`  
  Theme dictionary and helpers.
- `plot_theme_dotheatmap_pdf(...)` / `plot_theme_dotheatmap_from_results(...)`  
  Theme dot-heatmap.

## `tme_utils.R` — TME deconvolution helpers

- `undo_log_expr(expr, is_log, log_base)`  
  Reverse log transformation for TME tools that expect non-log input.
- `validate_tme_input(expr)`  
  Validate expression matrix: numeric, non-negative, rownames present; warn on Ensembl IDs.
- `deduplicate_expression_by_symbol(expr, symbol_col)`  
  Deduplicate rows by gene symbol, keeping the row with highest mean expression.
- `convert_mouse_symbols_to_human(expr, mouse_attr, human_attr, host, verbose)`  
  Convert mouse MGI symbols to human HGNC symbols via `biomaRt::getLDS`.
- `prepare_tme_expression(expr, is_log, species, log_base, verbose)`  
  One-stop preparation: undo log + optional mouse-to-human conversion.
- `calc_tme_barplot_size(n_samples, ...)` / `calc_tme_boxplot_size(n_celltypes, ...)`  
  Automatic figure-size calculation based on sample/cell-type counts.
- `run_iobr_deconvolution(expr, methods, perm, arrays, id_column)`  
  Run IOBR methods (`cibersort`, `epic`, `xcell`, `estimate`).
- `combine_tme_results(result_list, id_column)`  
  Inner-join IOBR result tables by sample ID.
- `melt_tme_results(tme_df, id_column, group_df, sample_col, group_col)`  
  Convert wide IOBR result to long format for plotting.
- `melt_estimate_scores(estimate_df, id_column, group_df, sample_col, group_col)`  
  Convert ESTIMATE scores to long format.
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

## `tme_utils.R` — TME deconvolution

- `undo_log_expr(expr, is_log, log_base)`  
  Reverse log transform.
- `run_iobr_deconvolution(expr, methods, ...)`  
  Run IOBR deconvolution methods.
- `combine_tme_results(result_list, id_column)`  
  Inner-join results across methods.
- `melt_tme_results(tme_df, ...)`  
  Long format for ggplot.
- `plot_tme_barplot_pdf(...)`  
  Stacked cell-fraction barplot.
- `plot_tme_boxplot_pdf(...)`  
  Cell-fraction box/violin plots.
- `melt_estimate_scores(estimate_df, ...)`  
  Long-format ESTIMATE scores.
- `plot_estimate_boxplot_pdf(...)`  
  ESTIMATE score boxplots.
