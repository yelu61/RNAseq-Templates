#!/usr/bin/env Rscript
# =============================================================================
# RNAseq_General — non-interactive analysis pipeline
# =============================================================================
# Runs the SAME standard workflow as notebooks/RNAseq_General.ipynb, driven by
# a config file instead of a parameter cell. Intended use:
#
#   1. Copy templates/General/ (config.R + run_analysis.R) into your project.
#   2. Edit config.R (paths, samples, groups, comparisons, thresholds, gene sets).
#   3. Run:  Rscript run_analysis.R
#      or:   Rscript run_analysis.R path/to/other_config.R
#
# The notebook remains the place for interactive exploration (adjusting gene
# sets for GSVA, ad-hoc visualizations). This script reproduces the standard
# run end-to-end and is what you would schedule / batch / re-run.
#
# LIB_DIR resolution: if ./RNAseq_lib exists next to this script it is used;
# otherwise the repository root is located via rprojroot and its RNAseq_lib is
# used. Set the environment variable RNASEQ_LIB_DIR to override explicitly.
# =============================================================================

options(stringsAsFactors = FALSE)

# ---- Locate this script and the config file ----------------------------------
invocation_dir <- normalizePath(getwd(), mustWork = TRUE)
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", cmd_args, value = TRUE))
script_dir <- if (length(file_arg) > 0 && nzchar(file_arg)) dirname(normalizePath(file_arg[1])) else invocation_dir

user_args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(user_args) >= 1) {
  normalizePath(path.expand(user_args[1]), mustWork = FALSE)
} else {
  file.path(script_dir, "config.R")
}
if (!file.exists(config_path)) {
  stop("Config file not found: ", config_path,
       "\nProvide one as: Rscript run_analysis.R path/to/config.R")
}
config_path <- normalizePath(config_path, mustWork = TRUE)

cat("========================================\n")
cat("RNAseq_General — run_analysis.R\n")
cat("Working dir :", getwd(), "\n")
cat("Config file :", config_path, "\n")
cat("========================================\n\n")

# Source the config into the global environment so all parameters are available.
source(config_path, local = globalenv())

# ---- Resolve RNAseq_lib -------------------------------------------------------
lib_dir <- Sys.getenv("RNASEQ_LIB_DIR", unset = NA_character_)
if (is.na(lib_dir) || !dir.exists(lib_dir)) {
  if (dir.exists(file.path(invocation_dir, "RNAseq_lib"))) {
    lib_dir <- file.path(invocation_dir, "RNAseq_lib")
  } else if (dir.exists(file.path(script_dir, "RNAseq_lib"))) {
    lib_dir <- file.path(script_dir, "RNAseq_lib")
  } else {
    repo_root <- tryCatch(rprojroot::find_root(rprojroot::is_git_root, path = script_dir), error = function(e) NA_character_)
    if (!is.na(repo_root) && dir.exists(file.path(repo_root, "RNAseq_lib"))) {
      lib_dir <- file.path(repo_root, "RNAseq_lib")
    } else if (dir.exists(file.path(dirname(script_dir), "RNAseq_lib"))) {
      lib_dir <- file.path(dirname(script_dir), "RNAseq_lib")
    } else {
      stop("Could not locate RNAseq_lib. Set RNASEQ_LIB_DIR or run from a project ",
           "that has RNAseq_lib/ alongside run_analysis.R.")
    }
  }
}
cat("RNAseq_lib  :", normalizePath(lib_dir), "\n\n")

# ---- Load libraries -----------------------------------------------------------
suppressPackageStartupMessages({
  library(DESeq2)
  library(clusterProfiler)
  if (SPECIES == "human") { library(org.Hs.eg.db); org_db <- org.Hs.eg.db } else { library(org.Mm.eg.db); org_db <- org.Mm.eg.db }
  library(ComplexHeatmap)
  library(circlize)
  library(matrixStats)
  library(EnhancedVolcano)
  library(GSVA)
  library(msigdbr)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(ggrepel)
  library(DOSE)
  library(enrichplot)
  library(ggpubr)
  library(tidyverse)
  library(data.table)
  library(pheatmap)
  library(ashr)
  library(BiocParallel)
  if (RUN_TF_ANALYSIS) {
    library(dorothea)
    library(viper)
    library(limma)
  }
})

# TME deconvolution is optional and has heavy/online dependencies. Source its
# helpers only when requested so a base run does not require IOBR/estimate.
if (isTRUE(RUN_TME)) {
  source(file.path(lib_dir, "tme_utils.R"))
  if (!exists("TME_ORTHOLOG_CACHE")) TME_ORTHOLOG_CACHE <- NULL  # configs predating this option
  if (isTRUE(RUN_TME_IOBR) && !requireNamespace("IOBR", quietly = TRUE)) {
    message("RUN_TME_IOBR is TRUE but the 'IOBR' package is not installed; skipping IOBR deconvolution.")
    RUN_TME_IOBR <- FALSE
  }
  if (isTRUE(RUN_TME_ESTIMATE) && !requireNamespace("estimate", quietly = TRUE)) {
    message("RUN_TME_ESTIMATE is TRUE but the 'estimate' package is not installed; skipping ESTIMATE.")
    RUN_TME_ESTIMATE <- FALSE
  }
}

for (f in c("plot_utils.R", "io_utils.R", "data_utils.R", "deg_utils.R",
            "enrichment_utils.R", "batch_utils.R", "design_utils.R", "report_utils.R")) {
  source(file.path(lib_dir, f))
}
# TME helpers are sourced later, only when RUN_TME is TRUE.

n_groups <- length(unique(GROUPS))
group_colors <- make_group_colors(GROUP_LEVELS)
colors_direction <- c("UP" = "#d6604d", "DOWN" = "#4393c3", "Not_Sig" = "#999999")
theme_set(theme_publication())

# ---- Output directories + config snapshot ------------------------------------
for (d in c("0-Config", "1-DEG", "2-GSEA", "3-Visualization")) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

config_objects <- c(
  "SPECIES", "INPUT_FILE", "INPUT_FORMAT", "GENE_NAME_COL", "BIOTYPE_COL", "BIOTYPE_FILTER",
  "COUNT_COLS", "SAMPLE_NAMES", "GROUPS", "GROUP_LEVELS", "BATCH_VECTOR", "PAIR_ID",
  "SAMPLE_EXCLUDE", "MIN_LIBRARY_SIZE", "MIN_DETECTED_GENES", "MAX_ZERO_FRACTION",
  "MIN_MEDIAN_CORRELATION_Z", "COMPARISONS", "DEG_PVALUE_COLUMN", "DEG_LFC_COLUMN",
  "GSEA_RANK_COLUMN", "THRESHOLD_GRID", "DEFAULT_THRESHOLD", "MIN_COUNT", "DESIGN_FORMULA",
  "PAIRWISE_TEST_METHOD", "PAIRWISE_P_ADJUST_METHOD", "RUN_TF_ANALYSIS", "RUN_COMPARECLUSTER",
  "COMPARECLUSTER_ONTOLOGY", "EXPORT_EXCEL", "GENERATE_HTML_REPORT",
  "RUN_TME", "TME_GENE_LENGTH_COLUMN", "TME_GENE_LENGTH_UNIT", "TME_GENE_START_COL",
  "TME_GENE_END_COL", "RUN_TME_ESTIMATE", "RUN_TME_IOBR", "TME_IOBR_METHODS",
  "TME_IOBR_PERM", "RUN_TME_SSGSEA", "TME_ORTHOLOG_CACHE"
)
config_lines <- c(
  "# RNAseq_General analysis configuration snapshot",
  paste0("# Saved: ", Sys.time()),
  unlist(lapply(config_objects, function(obj) {
    if (exists(obj, inherits = TRUE)) c(paste0("\n", obj, " <- "), capture.output(dput(get(obj, inherits = TRUE)))) else NULL
  }))
)
writeLines(config_lines, "./0-Config/analysis_config_used.R")

# =============================================================================
# 4. Data loading & preprocessing
# =============================================================================
rawcount <- read_count_table(INPUT_FILE, INPUT_FORMAT)
preprocess_summary_rows <- list()
preprocess_summary_rows[["raw_input"]] <- summarize_raw_input(rawcount, GENE_NAME_COL, BIOTYPE_COL, BIOTYPE_FILTER)

if (!is.null(BIOTYPE_COL) && BIOTYPE_COL %in% colnames(rawcount)) {
  rawcount <- rawcount[rawcount[[BIOTYPE_COL]] == BIOTYPE_FILTER, ]
}
preprocess_summary_rows[["after_biotype_filter"]] <- data.frame(step = "after_biotype_filter_rows", value = nrow(rawcount))

count_col_names <- detect_count_columns(rawcount, GENE_NAME_COL, COUNT_COLS)
validate_sample_design(SAMPLE_NAMES, GROUPS, GROUP_LEVELS, COMPARISONS, count_col_names)

countData <- build_count_matrix(rawcount, GENE_NAME_COL, count_col_names, SAMPLE_NAMES,
                                duplicate_report_file = "./1-DEG/Duplicated_gene_symbols.csv")
validate_count_matrix(countData)
countData_before_gene_filter <- countData
plot_count_distribution_pdf(countData_before_gene_filter, "./3-Visualization/Raw_count_distribution_boxplot.pdf")
preprocess_summary_rows[["after_gene_symbol_dedup"]] <- data.frame(step = "after_gene_symbol_dedup_rows", value = nrow(countData))

# colData, QC, optional exclusion, low-count filter
if (!is.null(PAIR_ID)) {
  validate_paired_design(SAMPLE_NAMES, GROUPS, PAIR_ID, GROUP_LEVELS)
  colData <- make_paired_col_data(SAMPLE_NAMES, GROUPS, GROUP_LEVELS, PAIR_ID)
  DESIGN_FORMULA <- build_paired_design_formula("condition", "pair_id")
  group <- colData$condition
  cat("Paired design enabled:", deparse(DESIGN_FORMULA), "\n")
} else {
  group <- factor(GROUPS, levels = GROUP_LEVELS)
  colData <- make_col_data(countData, SAMPLE_NAMES, GROUPS, GROUP_LEVELS)
}

# Batch covariate: surface diagnostics, and close the loop by writing batch into
# colData so a batch-aware DESIGN_FORMULA (e.g. ~ batch + condition) actually resolves.
if (!is.null(BATCH_VECTOR)) {
  validate_batch_design(BATCH_VECTOR, DESIGN_FORMULA, SAMPLE_NAMES)
  colData <- add_batch_col(colData, BATCH_VECTOR)
}

sample_qc <- calculate_sample_qc(countData, colData, group_col = "condition")
sample_qc <- flag_sample_qc(sample_qc,
                            min_library_size = MIN_LIBRARY_SIZE,
                            min_detected_genes = MIN_DETECTED_GENES,
                            max_zero_fraction = MAX_ZERO_FRACTION,
                            min_median_correlation_z = MIN_MEDIAN_CORRELATION_Z)
write.csv(sample_qc, "./3-Visualization/Sample_QC_metrics.csv", row.names = FALSE)
plot_sample_qc_pdf(sample_qc, "./3-Visualization/Sample_QC_metrics.pdf", group_colors)
plot_sample_correlation_pdf(sample_qc, "./3-Visualization/Sample_correlation_heatmap_raw_counts.pdf")

if (any(sample_qc$qc_flag)) {
  cat("QC-flagged samples (review before excluding):\n")
  print(sample_qc[sample_qc$qc_flag, c("sample", "qc_reason")])
}

if (length(SAMPLE_EXCLUDE) > 0) {
  exclusion_res <- apply_sample_exclusion(countData, colData, SAMPLE_EXCLUDE)
  countData <- exclusion_res$count_data
  colData <- exclusion_res$col_data
  group <- colData$condition
  cat("Excluded samples:", paste(exclusion_res$excluded_samples, collapse = ", "), "\n")
}

filter_res <- filter_low_count_genes(countData, group, MIN_COUNT)
countData <- filter_res$count_data
min_replicates <- filter_res$min_replicates

retention_df <- calculate_filter_retention(countData_before_gene_filter, countData)
write.csv(retention_df, "./3-Visualization/Filter_retention_by_sample.csv", row.names = FALSE)
plot_filter_retention_pdf(retention_df, "./3-Visualization/Filter_retention_by_sample.pdf")

preprocess_summary_rows[["after_sample_exclusion"]] <- data.frame(step = "samples_after_optional_exclusion", value = ncol(countData))
preprocess_summary_rows[["after_low_count_filter"]] <- data.frame(step = "genes_after_low_count_filter", value = nrow(countData))
preprocessing_summary <- write_preprocessing_summary("./1-DEG/Preprocessing_summary.csv", preprocess_summary_rows)

cat("After QC/exclusion/filtering:", nrow(countData), "genes x", ncol(countData), "samples\n")

# =============================================================================
# 5. DESeq2 object & QC
# =============================================================================
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = DESIGN_FORMULA)
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

plot_pca_pdf(vsd, GROUP_LEVELS, group_colors, "./3-Visualization/PCA_plot.pdf")

if (!is.null(BATCH_VECTOR) && length(BATCH_VECTOR) == ncol(vsd)) {
  plot_pca_by_batch_pdf(vsd, BATCH_VECTOR, "./3-Visualization/PCA_by_batch.pdf")
  pve_df <- summarize_pve_by_batch(vsd, BATCH_VECTOR, condition_vec = as.character(colData$condition))
  write.csv(pve_df, "./3-Visualization/Batch_PVE_summary.csv", row.names = FALSE)
  plot_batch_pve_pdf(pve_df, "./3-Visualization/Batch_PVE_barplot.pdf")
  cat("Batch-effect diagnostics saved.\n")
}

plot_sample_distance_pdf(vsd, "./3-Visualization/Sample_distance_heatmap.pdf")

write.csv(data.frame(gene_name = rownames(vsd_mat), vsd_mat, check.names = FALSE), "./1-DEG/vsd_matrix.csv", row.names = FALSE)
write.csv(as.data.frame(colData, check.names = FALSE), "./1-DEG/colData.csv", row.names = FALSE)
save(rawcount, countData_before_gene_filter, countData, group, colData, sample_qc,
     retention_df, preprocessing_summary, dds, vsd, vsd_mat, file = "./1-DEG/step0-input.Rdata")

# =============================================================================
# 6. DESeq2 differential expression
# =============================================================================
cat("Running DESeq2...\n")
dds <- DESeq(dds)

res_list <- extract_deseq2_results(dds = dds, comparisons = COMPARISONS, condition_col = "condition",
                                   alpha = max(THRESHOLD_GRID$p_cutoff), shrink_type = "ashr")
all_gene_deg_dir <- write_all_gene_deg_results(res_list, outdir = "./1-DEG")
deg_by_threshold <- build_deg_threshold_sets(res_list, THRESHOLD_GRID,
                                             pvalue_column = DEG_PVALUE_COLUMN, lfc_column = DEG_LFC_COLUMN)
default_res_list <- get_threshold_result(deg_by_threshold, DEFAULT_THRESHOLD)
deg_summary <- write_deg_threshold_outputs(deg_by_threshold, outdir = "./1-DEG",
                                           threshold_grid = THRESHOLD_GRID,
                                           pvalue_column = DEG_PVALUE_COLUMN, lfc_column = DEG_LFC_COLUMN)
lfc_strategy_summary <- write_lfc_strategy_summary(res_list, THRESHOLD_GRID, outdir = "./1-DEG",
                                                   pvalue_column = DEG_PVALUE_COLUMN)
deg_diagnostic_summary <- write_deg_diagnostic_summary(res_list, THRESHOLD_GRID, outdir = "./1-DEG")

if (isTRUE(EXPORT_EXCEL)) {
  write_deg_excel(deg_by_threshold, outdir = "./1-DEG", filename = "DEG_results.xlsx")
  cat("Excel export saved: ./1-DEG/DEG_results.xlsx\n")
}
save(res_list, deg_by_threshold, deg_summary, lfc_strategy_summary, deg_diagnostic_summary,
     all_gene_deg_dir, dds, vsd, vsd_mat, colData, GROUP_LEVELS, COMPARISONS, SPECIES,
     file = "./1-DEG/DEG_results.Rdata")
print(deg_summary)

# =============================================================================
# 7. DEG statistics visualization
# =============================================================================
plot_deg_summary_pdf(deg_summary, "./3-Visualization/DEG_threshold_summary_barplot.pdf")
if (exists("lfc_strategy_summary") && nrow(lfc_strategy_summary) > 0) {
  plot_lfc_strategy_summary_pdf(lfc_strategy_summary, "./3-Visualization/DEG_raw_vs_shrunken_summary_barplot.pdf")
}

# =============================================================================
# 8. Volcano plots (default threshold)
# =============================================================================
for (comp_name in names(res_list)) {
  plot_volcano_pdf(res_list[[comp_name]],
                   comp_name = paste0(comp_name, " (", DEFAULT_THRESHOLD, ")"),
                   pvalue_thresh = DEG_P_CUTOFF, log2fc_thresh = DEG_LFC_CUTOFF,
                   pvalue_column = DEG_PVALUE_COLUMN, lfc_column = DEG_LFC_COLUMN,
                   filename = paste0("./3-Visualization/Volcano_", comp_name, "_", DEFAULT_THRESHOLD, ".pdf"))
}

# =============================================================================
# 9. Heatmaps
# =============================================================================
top_genes_mad <- head(order(rowMads(vsd_mat), decreasing = TRUE), 1000)
plot_expression_heatmap_pdf(vsd_mat[top_genes_mad, ],
                            filename = "./3-Visualization/Heatmap_top1000_variable_genes.pdf",
                            title = "Top 1000 Variable Genes",
                            group = as.character(colData[colnames(vsd_mat), "condition"]),
                            group_levels = GROUP_LEVELS, group_colors = group_colors,
                            show_row_names = FALSE, show_column_names = FALSE,
                            cluster_columns = TRUE, width = mm_to_in(183), height = mm_to_in(247))

# Top DEGs heatmap
all_sig_genes <- unique(unlist(lapply(default_res_list, function(res) {
  res %>%
    filter(.data[[DEG_PVALUE_COLUMN]] < DEG_P_CUTOFF, abs(.data[[DEG_LFC_COLUMN]]) > DEG_LFC_CUTOFF) %>%
    arrange(.data[[DEG_PVALUE_COLUMN]]) %>%
    head(30) %>%
    pull(gene_name)
})))
sig_genes_in_mat <- intersect(all_sig_genes, rownames(vsd_mat))
if (length(sig_genes_in_mat) > 0) {
  plot_expression_heatmap_pdf(vsd_mat[sig_genes_in_mat, ],
                              filename = "./3-Visualization/Heatmap_topDEGs.pdf",
                              title = "Top Differentially Expressed Genes",
                              group = as.character(colData[colnames(vsd_mat), "condition"]),
                              group_levels = GROUP_LEVELS, group_colors = group_colors,
                              show_row_names = TRUE, show_column_names = TRUE,
                              row_font_size = 7, column_font_size = 8, cluster_columns = FALSE,
                              width = 10, height = max(8, 0.22 * length(sig_genes_in_mat) + 3))
}

# Key-gene + custom gene-set heatmaps
safe_plot_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)
key_genes_in_mat <- intersect(KEY_GENES, rownames(vsd_mat))
if (length(key_genes_in_mat) > 0) {
  plot_expression_heatmap_pdf(vsd_mat[key_genes_in_mat, , drop = FALSE],
                              filename = "./3-Visualization/Heatmap_key_genes.pdf",
                              title = "Key Genes",
                              group = as.character(colData[colnames(vsd_mat), "condition"]),
                              group_levels = GROUP_LEVELS, group_colors = group_colors,
                              show_row_names = TRUE, show_column_names = TRUE,
                              cluster_rows = TRUE, cluster_columns = TRUE, row_font_size = 8,
                              width = 8, height = max(5, 0.28 * length(key_genes_in_mat) + 3))
}

if (!is.null(custom_gene_sets) && length(custom_gene_sets) > 0) {
  dir.create("./3-Visualization/CustomGeneSets", showWarnings = FALSE, recursive = TRUE)
  row_upper <- toupper(rownames(vsd_mat))
  upper_to_real <- setNames(rownames(vsd_mat), row_upper)
  custom_heatmap_summary <- lapply(names(custom_gene_sets), function(set_name) {
    found <- unique(upper_to_real[intersect(toupper(custom_gene_sets[[set_name]]), names(upper_to_real))])
    data.frame(Signature = set_name, InputGenes = length(unique(custom_gene_sets[[set_name]])), FoundGenes = length(found))
  }) %>% bind_rows()
  write.csv(custom_heatmap_summary, "./3-Visualization/CustomGeneSets/custom_gene_set_heatmap_summary.csv", row.names = FALSE)
  for (set_name in names(custom_gene_sets)) {
    found <- unique(upper_to_real[intersect(toupper(custom_gene_sets[[set_name]]), names(upper_to_real))])
    if (length(found) < 2) next
    plot_expression_heatmap_pdf(vsd_mat[found, , drop = FALSE],
                                filename = file.path("./3-Visualization/CustomGeneSets", paste0("Heatmap_", safe_plot_name(set_name), ".pdf")),
                                title = paste0(set_name, " (", length(found), " genes)"),
                                group = as.character(colData[colnames(vsd_mat), "condition"]),
                                group_levels = GROUP_LEVELS, group_colors = group_colors,
                                show_row_names = TRUE, show_column_names = TRUE,
                                cluster_rows = TRUE, cluster_columns = TRUE, row_font_size = 8,
                                width = 8, height = max(4.5, 0.28 * length(found) + 3))
  }
}

# =============================================================================
# 10. Pathway enrichment (ORA + GSEA)
# =============================================================================
ranked_lists <- lapply(res_list, ranked_gene_list, rank_column = GSEA_RANK_COLUMN)
background_entrez <- map_symbols_to_entrez(rownames(countData), org_db)
background_universe <- unique(background_entrez$ENTREZID)
org_code <- ifelse(SPECIES == "human", "hsa", "mmu")

ora_threshold <- run_threshold_ora(res_list = res_list, threshold_grid = THRESHOLD_GRID,
                                   org_db = org_db, universe = background_universe, organism = org_code,
                                   pvalue_column = DEG_PVALUE_COLUMN, lfc_column = DEG_LFC_COLUMN,
                                   outdir = "./2-GSEA", plotdir = "./3-Visualization")
ora_summary <- ora_threshold$summary
write.csv(ora_summary, "./2-GSEA/ORA_threshold_summary.csv", row.names = FALSE)

# ORA summary plots
if (exists("ora_summary") && nrow(ora_summary) > 0) {
  ora_long <- ora_summary %>%
    pivot_longer(cols = c("GO_terms", "KEGG_terms"), names_to = "Database", values_to = "Terms")
  p_ora_summary <- ggplot(ora_long, aes(x = Comparison, y = Terms, fill = Database)) +
    geom_col(position = "dodge", width = 0.7) +
    facet_wrap(~ Threshold, scales = "free_x") +
    labs(x = NULL, y = "Significant ORA terms", title = "ORA Term Counts Across DEG Thresholds") +
    theme_publication(base_size = 8) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "top")
  save_pdf_plot(p_ora_summary, "./3-Visualization/ORA_threshold_summary_barplot.pdf", width = 10, height = 5)

  directional_cols <- intersect(c("GO_up_terms", "GO_down_terms", "KEGG_up_terms", "KEGG_down_terms"), colnames(ora_summary))
  if (length(directional_cols) > 0) {
    ora_direction_long <- ora_summary %>%
      pivot_longer(cols = all_of(directional_cols), names_to = "Category", values_to = "Terms") %>%
      separate(Category, into = c("Database", "Direction", "metric"), sep = "_", remove = FALSE)
    p_ora_direction <- ggplot(ora_direction_long, aes(x = Comparison, y = Terms, fill = Direction)) +
      geom_col(position = "dodge", width = 0.7) +
      facet_grid(Database ~ Threshold, scales = "free_x") +
      scale_fill_manual(values = c("up" = "#d6604d", "down" = "#4393c3")) +
      labs(x = NULL, y = "Directional ORA terms", title = "UP/DOWN ORA Term Counts") +
      theme_publication(base_size = 8) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "top")
    save_pdf_plot(p_ora_direction, "./3-Visualization/ORA_directional_summary_barplot.pdf", width = 11, height = 7)
  }
}

# GSEA (full ranked list)
gsea_results <- list()
gsea_root <- "./3-Visualization/GSEA"
gsea_overview_outdir <- file.path(gsea_root, "overview")
dir.create(gsea_overview_outdir, showWarnings = FALSE, recursive = TRUE)
for (comp_name in names(ranked_lists)) {
  entrezList <- make_entrez_ranked_list(ranked_lists[[comp_name]], org_db)
  if (length(entrezList) > 10) {
    ggo <- run_go_gsea(entrezList, org_db = org_db, ont = "ALL", min_size = 10, max_size = 500, p_cutoff = 0.05)
    if (!is.null(ggo) && nrow(as.data.frame(ggo)) > 0) {
      write.csv(as.data.frame(ggo), paste0("./2-GSEA/GSEA_GO_", comp_name, ".csv"), row.names = FALSE)
      plot_gsea_suite_pdf(ggo, file.path(gsea_overview_outdir, paste0("GSEA_GO_", comp_name)), paste("GSEA GO -", comp_name))
    }
    gkegg <- run_kegg_gsea(entrezList, organism = org_code, min_size = 10, max_size = 500, p_cutoff = 0.05)
    if (!is.null(gkegg) && nrow(as.data.frame(gkegg)) > 0) {
      write.csv(as.data.frame(gkegg), paste0("./2-GSEA/GSEA_KEGG_", comp_name, ".csv"), row.names = FALSE)
      plot_gsea_suite_pdf(gkegg, file.path(gsea_overview_outdir, paste0("GSEA_KEGG_", comp_name)), paste("GSEA KEGG -", comp_name))
    }
    gsea_results[[comp_name]] <- list(go = ggo, kegg = gkegg)
  }
}
# Cache the gseaResult objects so downstream/custom figures (e.g. single-term
# gseaplot2 in visualize_results.R) can be regenerated without re-running GSEA.
saveRDS(gsea_results, "./2-GSEA/gsea_results.rds")

# ---- Theme dot-heatmaps ----
theme_outdir <- file.path(gsea_root, "theme_maps")
dir.create(theme_outdir, showWarnings = FALSE, recursive = TRUE)
default_ora <- ora_threshold$results[[DEFAULT_THRESHOLD]]
go_ora_map <- lapply(default_ora, function(x) x$go)
kegg_ora_map <- lapply(default_ora, function(x) x$kegg)
go_gsea_map <- lapply(gsea_results, function(x) x$go)
kegg_gsea_map <- lapply(gsea_results, function(x) x$kegg)
drop_empty <- function(lst) lst[sapply(lst, function(x) !is.null(x) && nrow(as.data.frame(x)) > 0)]
go_ora_map <- drop_empty(go_ora_map); kegg_ora_map <- drop_empty(kegg_ora_map)
go_gsea_map <- drop_empty(go_gsea_map); kegg_gsea_map <- drop_empty(kegg_gsea_map)
theme_defs <- default_enrichment_themes()

if (length(go_ora_map) >= 1) {
  plot_theme_dotheatmap_from_results(go_ora_map, filename = file.path(theme_outdir, "Theme_dotheatmap_GO_ORA.pdf"),
                                     title = "GO ORA Biological Themes",
                                     subtitle = paste("GO-BP ORA |", DEFAULT_THRESHOLD, "| top terms per theme"),
                                     theme_defs = theme_defs, ontology_filter = "BP", top_n = 6)
}
if (length(kegg_ora_map) >= 1) {
  plot_theme_dotheatmap_from_results(kegg_ora_map, filename = file.path(theme_outdir, "Theme_dotheatmap_KEGG_ORA.pdf"),
                                     title = "KEGG ORA Pathway Themes",
                                     subtitle = paste("KEGG ORA |", DEFAULT_THRESHOLD, "| top terms per theme"),
                                     theme_defs = theme_defs, ontology_filter = NULL, top_n = 6)
}
if (length(go_gsea_map) >= 1) {
  plot_theme_dotheatmap_from_results(go_gsea_map, filename = file.path(theme_outdir, "Theme_dotheatmap_GO_GSEA.pdf"),
                                     title = "GO GSEA Biological Themes",
                                     subtitle = "GO-BP GSEA | top terms per theme",
                                     theme_defs = theme_defs, ontology_filter = "BP", top_n = 6)
}
if (length(kegg_gsea_map) >= 1) {
  plot_theme_dotheatmap_from_results(kegg_gsea_map, filename = file.path(theme_outdir, "Theme_dotheatmap_KEGG_GSEA.pdf"),
                                     title = "KEGG GSEA Pathway Themes",
                                     subtitle = "KEGG GSEA | top terms per theme",
                                     theme_defs = theme_defs, ontology_filter = NULL, top_n = 6)
}

# ---- Single-term GSEA figures ----
single_term_outdir <- file.path(gsea_root, "running_curves")
dir.create(single_term_outdir, showWarnings = FALSE, recursive = TRUE)
for (comp_name in names(gsea_results)) {
  ggo <- gsea_results[[comp_name]]$go
  gkegg <- gsea_results[[comp_name]]$kegg
  if (!is.null(ggo) && nrow(as.data.frame(ggo)) > 0) {
    ggo_df <- as.data.frame(ggo); ggo_df <- ggo_df[order(ggo_df$p.adjust, -abs(ggo_df$NES)), ]
    top_terms <- rbind(utils::head(ggo_df[ggo_df$NES > 0, ], 3), utils::head(ggo_df[ggo_df$NES < 0, ], 3))
    plot_gsea_term_figures_from_df(ggo, top_terms, outdir = file.path(single_term_outdir, paste0("GO_", comp_name)),
                                   contrast_label = comp_name, prefix = "gseaplot2_GO")
  }
  if (!is.null(gkegg) && nrow(as.data.frame(gkegg)) > 0) {
    gkegg_df <- as.data.frame(gkegg); gkegg_df <- gkegg_df[order(gkegg_df$p.adjust, -abs(gkegg_df$NES)), ]
    top_terms <- rbind(utils::head(gkegg_df[gkegg_df$NES > 0, ], 3), utils::head(gkegg_df[gkegg_df$NES < 0, ], 3))
    plot_gsea_term_figures_from_df(gkegg, top_terms, outdir = file.path(single_term_outdir, paste0("KEGG_", comp_name)),
                                   contrast_label = comp_name, prefix = "gseaplot2_KEGG")
  }
}

# ---- compareCluster (multi-group) ----
if (RUN_COMPARECLUSTER && n_groups >= 3) {
  cp_up <- list(); cp_down <- list()
  for (comp_name in names(res_list)) {
    ed <- genes_for_enrichment(res_list[[comp_name]], DEG_P_CUTOFF, DEG_LFC_CUTOFF,
                               pvalue_column = DEG_PVALUE_COLUMN, lfc_column = DEG_LFC_COLUMN)
    if (length(ed$up) >= 5) {
      up_entrez <- map_symbols_to_entrez(ed$up, org_db)
      if (nrow(up_entrez) >= 5) cp_up[[comp_name]] <- up_entrez$ENTREZID
    }
    if (length(ed$down) >= 5) {
      down_entrez <- map_symbols_to_entrez(ed$down, org_db)
      if (nrow(down_entrez) >= 5) cp_down[[comp_name]] <- down_entrez$ENTREZID
    }
  }

  run_comparecluster_go <- function(gene_lists, direction_label) {
    if (length(gene_lists) < 2) return(NULL)
    mg_go <- compareCluster(gene_lists, fun = "enrichGO", OrgDb = org_db, universe = background_universe,
                            ont = COMPARECLUSTER_ONTOLOGY, pAdjustMethod = "BH", pvalueCutoff = 0.1, readable = TRUE)
    if (!is.null(mg_go) && nrow(as.data.frame(mg_go)) > 0) {
      write.csv(as.data.frame(mg_go), paste0("./2-GSEA/CompareCluster_GO_", direction_label, ".csv"), row.names = FALSE)
      plot_comparecluster_dotplot_pdf(mg_go, paste0("./3-Visualization/CompareCluster_GO_", direction_label, ".pdf"),
                                      paste("GO Enrichment -", direction_label, "genes"), show_category = 10, width = 10, height = 12)
    }
    mg_go
  }
  run_comparecluster_kegg <- function(gene_lists, direction_label) {
    if (length(gene_lists) < 2) return(NULL)
    mg_kegg <- compareCluster(gene_lists, fun = "enrichKEGG", organism = org_code, universe = background_universe,
                              pAdjustMethod = "BH", pvalueCutoff = 0.1)
    if (!is.null(mg_kegg) && nrow(as.data.frame(mg_kegg)) > 0) {
      write.csv(as.data.frame(mg_kegg), paste0("./2-GSEA/CompareCluster_KEGG_", direction_label, ".csv"), row.names = FALSE)
      plot_comparecluster_dotplot_pdf(mg_kegg, paste0("./3-Visualization/CompareCluster_KEGG_", direction_label, ".pdf"),
                                      paste("KEGG Enrichment -", direction_label, "genes"), show_category = 10, width = 10, height = 10)
    }
    mg_kegg
  }
  run_comparecluster_go(cp_up, "up"); run_comparecluster_go(cp_down, "down")
  run_comparecluster_kegg(cp_up, "up"); run_comparecluster_kegg(cp_down, "down")
} else if (RUN_COMPARECLUSTER) {
  cat("CompareCluster skipped: only", n_groups, "group(s) < 3.\n")
}

# =============================================================================
# 11. GSVA
# =============================================================================
gs_for_gsva <- list()
gsva_scores_df <- NULL
if (!is.null(custom_gene_sets)) {
  row_upper <- toupper(rownames(vsd_mat))
  upper_to_real <- setNames(rownames(vsd_mat), row_upper)
  for (gs_name in names(custom_gene_sets)) {
    found <- unique(upper_to_real[intersect(toupper(custom_gene_sets[[gs_name]]), names(upper_to_real))])
    if (length(found) >= 3) gs_for_gsva[[gs_name]] <- found
  }
  if (length(gs_for_gsva) > 0) {
    params <- gsvaParam(as.matrix(vsd_mat), gs_for_gsva, minSize = 1, maxSize = Inf,
                        kcdf = "Gaussian", tau = 1, maxDiff = TRUE)
    cat("Running GSVA...\n")
    gsva_scores <- gsva(params, verbose = FALSE, BPPARAM = SerialParam())
    gsva_scores_df <- as.data.frame(gsva_scores)
    write.csv(gsva_scores_df, "./2-GSEA/GSVA_scores.csv")

    plot_expression_heatmap_pdf(gsva_scores_df, filename = "./3-Visualization/GSVA_heatmap.pdf",
                                title = "GSVA Scores",
                                group = as.character(colData[colnames(gsva_scores_df), "condition"]),
                                group_levels = GROUP_LEVELS, group_colors = group_colors,
                                scale_rows = FALSE, show_row_names = TRUE, show_column_names = FALSE,
                                row_font_size = 10, cluster_rows = TRUE, cluster_columns = TRUE,
                                width = 8, height = max(3, length(gs_for_gsva) * 0.8), heatmap_name = "GSVA Score")

    sample_conditions <- data.frame(
      sample = colnames(gsva_scores_df),
      condition = factor(as.character(colData[colnames(gsva_scores_df), "condition"]), levels = GROUP_LEVELS))
    gsva_long <- gsva_scores_df %>%
      as.data.frame() %>%
      mutate(Signature = rownames(gsva_scores_df)) %>%
      pivot_longer(cols = -Signature, names_to = "sample", values_to = "Score") %>%
      left_join(sample_conditions, by = "sample")
    plot_group_boxplot_pdf(gsva_long, value_col = "Score", group_col = "condition", facet_col = "Signature",
                           comparisons = combn(GROUP_LEVELS, 2, simplify = FALSE),
                           method = PAIRWISE_TEST_METHOD, p_adjust_method = PAIRWISE_P_ADJUST_METHOD,
                           title = "GSVA Scores", ylab = "Score", group_colors = group_colors,
                           filename = "./3-Visualization/GSVA_boxplot.pdf", width = 12,
                           height = 3 * ceiling(length(gs_for_gsva) / 2))
    dir.create("./3-Visualization/GSVA_boxplots", showWarnings = FALSE, recursive = TRUE)
    for (sig in unique(gsva_long$Signature)) {
      sig_data <- gsva_long %>% filter(Signature == sig)
      plot_group_violin_boxplot_pdf(sig_data, value_col = "Score", group_col = "condition",
                                    comparisons = combn(GROUP_LEVELS, 2, simplify = FALSE),
                                    method = PAIRWISE_TEST_METHOD, p_adjust_method = PAIRWISE_P_ADJUST_METHOD,
                                    title = sig, ylab = "GSVA Score", group_colors = group_colors,
                                    filename = file.path("./3-Visualization/GSVA_boxplots", paste0("GSVA_boxplot_", safe_plot_name(sig), ".pdf")),
                                    width = 5.5, height = 6)
    }
  }
}

# =============================================================================
# 12. Single-gene expression plots
# =============================================================================
single_gene_dir <- "./3-Visualization/SingleGene"
dir.create(single_gene_dir, showWarnings = FALSE, recursive = TRUE)
plot_gene_expression <- function(gene_name, vsd_mat_obj, dds_obj) {
  if (!gene_name %in% rownames(vsd_mat_obj)) return(NULL)
  plot_data <- data.frame(
    sample = colnames(vsd_mat_obj),
    expression = as.numeric(vsd_mat_obj[gene_name, colnames(vsd_mat_obj)]),
    condition = as.character(colData(dds_obj)[colnames(vsd_mat_obj), "condition"]))
  plot_data$condition <- factor(plot_data$condition, levels = GROUP_LEVELS)
  plot_group_bar_sem_pdf(plot_data, value_col = "expression", group_col = "condition",
                         comparisons = combn(GROUP_LEVELS, 2, simplify = FALSE),
                         method = PAIRWISE_TEST_METHOD, p_adjust_method = PAIRWISE_P_ADJUST_METHOD,
                         title = gene_name, ylab = "VST Expression", group_colors = group_colors,
                         filename = file.path(single_gene_dir, paste0("DEG_barplot_", gene_name, ".pdf")),
                         width = 5.5, height = 6, y_from_zero = TRUE)
}
found_genes <- intersect(KEY_GENES, rownames(vsd_mat))
write.csv(data.frame(gene = found_genes), file.path(single_gene_dir, "key_genes_plotted.csv"), row.names = FALSE)
for (gene in found_genes) plot_gene_expression(gene, vsd_mat, dds)

# =============================================================================
# 12.5 DEG set overlap
# =============================================================================
RUN_DEG_OVERLAP <- length(COMPARISONS) >= 2 || nrow(THRESHOLD_GRID) > 1
if (RUN_DEG_OVERLAP) {
  # Overlap figures are diagnostic; a single failure must not abort the whole run
  # (and skip Analysis_summary.txt) after the expensive DEG/GSEA/GSVA steps.
  tryCatch({
    sig_sets <- build_deg_gene_sets(default_res_list, pvalue_thresh = DEG_P_CUTOFF, log2fc_thresh = DEG_LFC_CUTOFF,
                                    pvalue_column = DEG_PVALUE_COLUMN, lfc_column = DEG_LFC_COLUMN, direction = "sig")
    if (length(sig_sets) >= 2) {
      plot_deg_upset_pdf(sig_sets, filename = "./3-Visualization/DEG_overlap_upset.pdf",
                         title = paste("DEG Overlaps -", DEFAULT_THRESHOLD), width = 14, height = 8)
      plot_deg_overlap_heatmap_pdf(sig_sets, filename = "./3-Visualization/DEG_overlap_jaccard_heatmap.pdf",
                                   title = paste("DEG Set Overlap (Jaccard) -", DEFAULT_THRESHOLD))
      common_genes <- extract_intersection_genes(sig_sets, mode = "all")
      write.csv(data.frame(gene = common_genes), "./3-Visualization/DEG_overlap_common_genes.csv", row.names = FALSE)
    }
    up_sets <- build_deg_gene_sets(default_res_list, pvalue_thresh = DEG_P_CUTOFF, log2fc_thresh = DEG_LFC_CUTOFF,
                                   pvalue_column = DEG_PVALUE_COLUMN, lfc_column = DEG_LFC_COLUMN, direction = "up")
    down_sets <- build_deg_gene_sets(default_res_list, pvalue_thresh = DEG_P_CUTOFF, log2fc_thresh = DEG_LFC_CUTOFF,
                                     pvalue_column = DEG_PVALUE_COLUMN, lfc_column = DEG_LFC_COLUMN, direction = "down")
    if (length(up_sets) >= 2) plot_deg_upset_pdf(up_sets, filename = "./3-Visualization/DEG_UP_overlap_upset.pdf", title = "UP DEG Overlaps")
    if (length(down_sets) >= 2) plot_deg_upset_pdf(down_sets, filename = "./3-Visualization/DEG_DOWN_overlap_upset.pdf", title = "DOWN DEG Overlaps")
  }, error = function(e) message("Skipping DEG overlap figures: ", conditionMessage(e)))
}

# =============================================================================
# 13.5 Tumor microenvironment deconvolution (optional; needs TPM, not VST)
# =============================================================================
if (isTRUE(RUN_TME)) {
  cat("Running TME deconvolution (requires gene lengths for TPM)...\n")
  tme_outdir <- "4-TME"
  dir.create(tme_outdir, showWarnings = FALSE, recursive = TRUE)

  # Build TPM from raw counts + gene lengths. VST cannot be inverted to TPM.
  gene_lengths_kb <- extract_gene_lengths(
    rawcount,
    id_col = GENE_NAME_COL,
    length_col = TME_GENE_LENGTH_COLUMN,
    start_col = TME_GENE_START_COL,
    end_col = TME_GENE_END_COL,
    length_unit = TME_GENE_LENGTH_UNIT
  )
  expr_tpm <- counts_to_tpm(countData, gene_lengths_kb)
  expr_tpm <- expr_tpm[, colnames(countData), drop = FALSE]
  write.csv(data.frame(gene_name = rownames(expr_tpm), expr_tpm, check.names = FALSE),
            file.path(tme_outdir, "TPM_matrix.csv"), row.names = FALSE)

  # Prepare (dedup; mouse->human ortholog conversion when SPECIES == "mouse").
  # ortholog_cache makes repeat runs offline-deterministic (NULL = always online).
  expr_tme <- prepare_tme_expression(as.data.frame(expr_tpm, check.names = FALSE),
                                     is_log = FALSE, species = SPECIES,
                                     ortholog_cache = TME_ORTHOLOG_CACHE, verbose = TRUE)
  validate_tme_input(expr_tme)

  tme_meta <- data.frame(sample = colnames(countData),
                         condition = as.character(colData[colnames(countData), "condition"]),
                         check.names = FALSE)
  group_df_for_plot <- tme_meta[, c("sample", "condition"), drop = FALSE]

  # ---- Native ESTIMATE ----
  if (isTRUE(RUN_TME_ESTIMATE)) {
    utils::data("common_genes", package = "estimate", envir = environment())
    utils::data("SI_geneset", package = "estimate", envir = environment())
    estimate_df <- data.frame(NAME = rownames(expr_tme), Description = NA, expr_tme, check.names = FALSE)
    write.table(estimate_df, file.path(tme_outdir, "estimate_input.gct"), sep = "\t", quote = FALSE, row.names = FALSE)
    estimate::filterCommonGenes(input.f = file.path(tme_outdir, "estimate_input.gct"),
                                output.f = file.path(tme_outdir, "estimate_common_genes.gct"), id = "GeneSymbol")
    estimate::estimateScore(file.path(tme_outdir, "estimate_common_genes.gct"),
                            file.path(tme_outdir, "estimate_scores.gct"), platform = "illumina")
    estimate_scores <- read.table(file.path(tme_outdir, "estimate_scores.gct"), skip = 2, header = TRUE, sep = "\t", check.names = FALSE)
    rownames(estimate_scores) <- estimate_scores$NAME
    # Keep only the score rows and the real sample columns. ESTIMATE writes extra
    # Description/Description.1 columns that are not samples; select explicitly.
    score_rows <- setdiff(rownames(estimate_scores), c("NAME", "Description"))
    sample_cols <- intersect(colnames(estimate_scores), colnames(countData))
    estimate_scores <- as.data.frame(t(estimate_scores[score_rows, sample_cols, drop = FALSE]))
    estimate_scores$name <- rownames(estimate_scores)
    estimate_scores <- estimate_scores[, c("name", score_rows), drop = FALSE]
    write.csv(estimate_scores, file.path(tme_outdir, "ESTIMATE_scores.csv"), row.names = FALSE)
    cat("ESTIMATE scores saved.\n")
  }

  # ---- IOBR multi-algorithm deconvolution ----
  iobr_results <- list()
  if (isTRUE(RUN_TME_IOBR)) {
    iobr_results <- run_iobr_deconvolution(expr_tme, methods = TME_IOBR_METHODS,
                                           perm = TME_IOBR_PERM, arrays = FALSE, id_column = "sample")
    for (method in names(iobr_results)) {
      write.csv(iobr_results[[method]], file.path(tme_outdir, paste0("IOBR_", method, ".csv")), row.names = FALSE)
    }
    if (length(iobr_results) >= 2) {
      write.csv(combine_tme_results(iobr_results, id_column = "sample"),
                file.path(tme_outdir, "IOBR_TME_combined.csv"), row.names = FALSE)
    }

    if ("estimate" %in% names(iobr_results)) {
      est_long <- melt_estimate_scores(iobr_results[["estimate"]], id_column = "sample",
                                       group_df = group_df_for_plot, sample_col = "sample", group_col = "condition")
      plot_estimate_boxplot_pdf(est_long, group_col = "condition",
                                filename = file.path(tme_outdir, "IOBR_ESTIMATE_scores_boxplot.pdf"),
                                title = "ESTIMATE Scores by Group", group_colors = group_colors,
                                save_individual = TRUE, individual_prefix = file.path(tme_outdir, "IOBR_ESTIMATE"))
      plot_tme_heatmap_pdf(iobr_results[["estimate"]], tme_meta, group_col = "condition",
                           sample_col = "sample", group_colors = group_colors,
                           filename = file.path(tme_outdir, "IOBR_ESTIMATE_heatmap.pdf"), title = "IOBR ESTIMATE Scores")
    }

    for (method_name in intersect(c("cibersort", "epic"), names(iobr_results))) {
      method_label <- toupper(method_name)
      method_long <- melt_tme_results(iobr_results[[method_name]], id_column = "sample",
                                      group_df = group_df_for_plot, sample_col = "sample", group_col = "condition") |>
        dplyr::filter(!grepl("P-value|Correlation|RMSE", .data$cell_type))
      bar_size <- calc_tme_barplot_size(length(unique(method_long$sample)), length(unique(method_long$cell_type)))
      box_size <- calc_tme_boxplot_size(length(unique(method_long$cell_type)))
      prefix <- file.path(tme_outdir, paste0("IOBR_", method_label))
      plot_tme_barplot_pdf(method_long, group_col = "condition", sample_col = "sample",
                           filename = paste0(prefix, "_barplot.pdf"), title = paste(method_label, "Cell Fractions"),
                           width = bar_size["width"], height = bar_size["height"], group_colors = group_colors)
      plot_tme_boxplot_pdf(method_long, group_col = "condition", value_col = "fraction",
                           filename = paste0(prefix, "_boxplot.pdf"), title = paste(method_label, "Cell Fractions by Group"),
                           width = box_size["width"], height = box_size["height"], group_colors = group_colors)
    }

    if ("xcell" %in% names(iobr_results)) {
      plot_tme_heatmap_pdf(iobr_results[["xcell"]], tme_meta, group_col = "condition",
                           sample_col = "sample", group_colors = group_colors,
                           filename = file.path(tme_outdir, "IOBR_xCell_heatmap.pdf"), title = "xCell Scores",
                           width = 10, height = 12)
    }
  }

  # ---- ssGSEA immune signatures ----
  if (isTRUE(RUN_TME_SSGSEA)) {
    expr_for_ssgsea <- as.data.frame(expr_tme, check.names = FALSE)
    row_upper <- toupper(rownames(expr_for_ssgsea))
    upper_to_real <- setNames(rownames(expr_for_ssgsea), row_upper)
    gs <- lapply(immune_gene_sets, function(x) unique(upper_to_real[intersect(toupper(x), names(upper_to_real))]))
    gs <- gs[lengths(gs) >= 2]
    if (length(gs) > 0) {
      params <- gsvaParam(as.matrix(expr_for_ssgsea), gs, kcdf = "Gaussian", minSize = 2, maxSize = Inf)
      ssgsea_scores <- gsva(params, verbose = FALSE)
      write.csv(ssgsea_scores, file.path(tme_outdir, "ssGSEA_immune_scores.csv"))

      score_df <- as.data.frame(t(ssgsea_scores)) %>% rownames_to_column("sample") %>%
        left_join(tme_meta, by = "sample")
      score_long <- score_df %>% pivot_longer(cols = names(gs), names_to = "signature", values_to = "score")
      score_long$condition <- factor(score_long$condition, levels = GROUP_LEVELS)
      plot_group_boxplot_pdf(score_long, value_col = "score", group_col = "condition", facet_col = "signature",
                             comparisons = if (length(GROUP_LEVELS) >= 2) combn(GROUP_LEVELS, 2, simplify = FALSE) else list(),
                             method = PAIRWISE_TEST_METHOD, p_adjust_method = PAIRWISE_P_ADJUST_METHOD,
                             title = "ssGSEA Immune Signature Scores", ylab = "ssGSEA score",
                             group_colors = group_colors,
                             filename = file.path(tme_outdir, "ssGSEA_group_boxplot.pdf"), width = 12, height = 8)
      cat("ssGSEA immune signature scores saved.\n")
    } else {
      cat("ssGSEA skipped: no immune signature had >= 2 matched genes.\n")
    }
  }

  cat("TME outputs written to", tme_outdir, "\n")
}

# =============================================================================
# 14. Transcription factor activity (optional)
# =============================================================================
if (RUN_TF_ANALYSIS) {
  if (SPECIES == "human") { data(dorothea_hs, package = "dorothea"); regulon_df <- dorothea_hs } else { data(dorothea_mm, package = "dorothea"); regulon_df <- dorothea_mm }
  regulon_df <- regulon_df %>% filter(confidence %in% c("A", "B", "C"))
  regulon_list <- dorothea::df2regulon(regulon_df)
  tf_activity_matrix <- viper(eset = vsd_mat, regulon = regulon_list, nes = TRUE, method = "none", verbose = FALSE)
  write.csv(tf_activity_matrix, "./2-GSEA/TF_activity_matrix.csv")

  tf_design <- model.matrix(~ 0 + vsd$condition)
  colnames(tf_design) <- levels(vsd$condition)
  tf_res_list <- list()
  for (comp in COMPARISONS) {
    comp_name <- comp[1]; treat <- comp[2]; ctrl <- comp[3]
    if (!(treat %in% colnames(tf_design)) || !(ctrl %in% colnames(tf_design))) next
    contrast_mat <- makeContrasts(contrasts = paste0(treat, " - ", ctrl), levels = tf_design)
    fit2 <- eBayes(contrasts.fit(lmFit(tf_activity_matrix, tf_design), contrast_mat))
    diff_tfs <- topTable(fit2, number = Inf, sort.by = "P")
    diff_tfs$TF <- rownames(diff_tfs)
    diff_tfs$significance <- ifelse(diff_tfs$adj.P.Val < PADJ_THRESH & abs(diff_tfs$logFC) > LOG2FC_THRESH,
                                    ifelse(diff_tfs$logFC > 0, "Up", "Down"), "Not_Sig")
    tf_res_list[[comp_name]] <- diff_tfs
    write.csv(diff_tfs, paste0("./2-GSEA/TF_diff_activity_", comp_name, ".csv"), row.names = FALSE)
  }

  if (length(tf_res_list) > 0) {
    top_tfs <- unique(unlist(lapply(tf_res_list, function(df) df %>% filter(adj.P.Val < PADJ_THRESH, abs(logFC) > LOG2FC_THRESH) %>% head(15) %>% pull(TF))))
    if (length(top_tfs) > 0) {
      tf_heatmap_mat <- tf_activity_matrix[intersect(top_tfs, rownames(tf_activity_matrix)), ]
      plot_expression_heatmap_pdf(tf_heatmap_mat, filename = "./3-Visualization/TF_activity_heatmap.pdf",
                                  title = "TF Activity", group = as.character(colData[colnames(tf_heatmap_mat), "condition"]),
                                  group_levels = GROUP_LEVELS, group_colors = group_colors, scale_rows = FALSE,
                                  show_row_names = TRUE, show_column_names = FALSE, row_font_size = 9,
                                  cluster_rows = TRUE, cluster_columns = TRUE, width = 8,
                                  height = max(4, length(top_tfs) * 0.3), heatmap_name = "TF Activity (NES)")
    }
    comp_name <- names(tf_res_list)[1]
    diff_tfs <- tf_res_list[[comp_name]]
    top_tf <- rbind(diff_tfs %>% filter(adj.P.Val < PADJ_THRESH, logFC > LOG2FC_THRESH) %>% head(10),
                    diff_tfs %>% filter(adj.P.Val < PADJ_THRESH, logFC < -LOG2FC_THRESH) %>% head(10))
    if (nrow(top_tf) > 0) {
      top_tf$Regulation <- ifelse(top_tf$logFC > 0, "Up", "Down")
      p_tf_bar <- ggplot(top_tf, aes(x = reorder(TF, logFC), y = logFC, fill = Regulation)) +
        geom_bar(stat = "identity", width = 0.7) + coord_flip() +
        scale_fill_manual(values = colors_direction) +
        labs(title = paste("Top Differentially Active TFs -", comp_name), x = NULL, y = "Log2 Fold Change (TF Activity)") +
        theme_publication(base_size = 8) + theme(legend.position = "top")
      save_pdf_plot(
        p_tf_bar,
        paste0("./3-Visualization/TF_barplot_", comp_name, ".pdf"),
        width = mm_to_in(183),
        height = min(mm_to_in(247), max(4.8, 0.34 * nrow(top_tf) + 2.2))
      )
    }
  }
}

# =============================================================================
# 13. Summary report (text)
# =============================================================================
summary_report <- paste0(
  "========================================\n",
  "RNA-seq Analysis Summary Report\n",
  "========================================\n\n",
  "1. Data Overview\n",
  "   - Species: ", SPECIES, "\n",
  "   - Samples after optional QC exclusion: ", ncol(countData), "\n",
  "   - Samples manually excluded: ", ifelse(length(SAMPLE_EXCLUDE) == 0, "None", paste(SAMPLE_EXCLUDE, collapse = ", ")), "\n",
  "   - Groups after optional QC exclusion: ", paste(unique(as.character(colData$condition)), collapse = ", "), "\n",
  "   - DESeq2 design: ", deparse(DESIGN_FORMULA), "\n",
  "   - DEG p-value column: ", DEG_PVALUE_COLUMN, "\n",
  "   - DEG log2FC column: ", DEG_LFC_COLUMN, "\n",
  "   - GSEA rank column: ", GSEA_RANK_COLUMN, "\n",
  "   - Pairwise plot test: ", PAIRWISE_TEST_METHOD, " with ", PAIRWISE_P_ADJUST_METHOD, " adjustment\n",
  "   - Low-count filter: count >= ", MIN_COUNT, " in at least ", min_replicates, " samples\n",
  "   - Genes after filtering: ", nrow(countData), "\n\n",
  "2. Differential Expression Default Threshold: ", DEFAULT_THRESHOLD,
  " (", DEG_PVALUE_COLUMN, " < ", DEG_P_CUTOFF, " & |", DEG_LFC_COLUMN, "| > ", DEG_LFC_CUTOFF, ")\n"
)
for (comp_name in names(default_res_list)) {
  res <- default_res_list[[comp_name]]
  n_up <- sum(res$significance == "Up", na.rm = TRUE)
  n_down <- sum(res$significance == "Down", na.rm = TRUE)
  summary_report <- paste0(summary_report, "   - ", comp_name, ": ", n_up + n_down,
                           " DEGs (Up: ", n_up, ", Down: ", n_down, ")\n")
}
summary_report <- paste0(summary_report, "\n3. Output Files\n",
                         "   - ./0-Config/analysis_config_used.R\n",
                         "   - ./1-DEG/ (all-gene, threshold, diagnostics, Excel, Rdata)\n",
                         "   - ./2-GSEA/ (ORA/GSEA/GSVA/TF results)\n",
                         "   - ./3-Visualization/GSEA/ (overview, running curves, theme maps)\n",
                         "   - ./3-Visualization/SingleGene/ (single-gene expression plots)\n",
                         "   - ./3-Visualization/ (QC, DEG, ORA and other PDF figures)\n")
if (isTRUE(RUN_TME)) {
  summary_report <- paste0(summary_report,
                           "   - ./4-TME/ (TPM matrix, ESTIMATE/IOBR/ssGSEA deconvolution)\n")
}
summary_report <- paste0(summary_report, "\n========================================\n")
cat(summary_report)
writeLines(summary_report, "./Analysis_summary.txt")
writeLines(capture.output(sessionInfo()), "./sessionInfo.txt")

# =============================================================================
# 15. HTML report
# =============================================================================
if (isTRUE(GENERATE_HTML_REPORT)) {
  has_quarto <- requireNamespace("quarto", quietly = TRUE) && nzchar(Sys.which("quarto"))
  has_rmarkdown <- requireNamespace("rmarkdown", quietly = TRUE)
  if (!has_quarto && !has_rmarkdown) {
    message("Skipping HTML report: install the quarto CLI (+ R 'quarto' package) or the 'rmarkdown' package.")
  } else {
    report_path <- tryCatch(
      render_analysis_report(outdir = ".", report_file = "RNAseq_report.html",
                             params = list(title = REPORT_TITLE, author = Sys.info()[["user"]])),
      error = function(e) { message("HTML report failed: ", conditionMessage(e)); NULL }
    )
    if (!is.null(report_path)) cat("HTML report written to:", report_path, "\n")
  }
}

cat("\n========================================\n")
cat("RNAseq_General analysis COMPLETE.\n")
cat("Outputs: 1-DEG/  2-GSEA/  3-Visualization/")
if (isTRUE(RUN_TME)) cat("  4-TME/")
cat("  Analysis_summary.txt")
if (isTRUE(GENERATE_HTML_REPORT)) cat("  RNAseq_report.html")
cat("\n========================================\n")
