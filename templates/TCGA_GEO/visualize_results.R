#!/usr/bin/env Rscript
# =============================================================================
# RNAseq_TCGA_GEO — targeted re-visualization from saved results
# =============================================================================
# Run this AFTER run_analysis.R has produced results in the current project.
# It reads the saved intermediate objects — no DESeq2 / ORA / GSEA / survival
# recompute:
#
#   1-DEG/TCGA_GEO_results.Rdata   res_df, expr_log, clinical, surv_df,
#                                  uni_cox, multi_cox, ego, ekegg, gsea_go, ...
#
# Use it to iterate cheaply on figures: restyle a volcano, re-draw a single-gene
# expression boxplot, rebuild KM / Cox forest plots, or a theme dot-heatmap —
# without re-running the pipeline.
#
#   Rscript visualize_results.R
#
# Edit the CONFIG block below. Each section is independent and guarded, so a
# missing input simply skips that section with a message.
# =============================================================================

options(stringsAsFactors = FALSE)

# ---- Resolve project dir + RNAseq_lib ---------------------------------------
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", cmd_args, value = TRUE))
script_dir <- if (length(file_arg) > 0 && nzchar(file_arg)) dirname(normalizePath(file_arg[1])) else getwd()
setwd(script_dir)

lib_dir <- Sys.getenv("RNASEQ_LIB_DIR", unset = NA_character_)
if (is.na(lib_dir) || !dir.exists(lib_dir)) {
  if (dir.exists("RNAseq_lib")) {
    lib_dir <- "RNAseq_lib"
  } else {
    repo_root <- tryCatch(rprojroot::find_root(rprojroot::is_git_root), error = function(e) NA_character_)
    if (!is.na(repo_root) && dir.exists(file.path(repo_root, "RNAseq_lib"))) lib_dir <- file.path(repo_root, "RNAseq_lib")
    else if (dir.exists(file.path("..", "RNAseq_lib"))) lib_dir <- file.path("..", "RNAseq_lib")
    else stop("Could not locate RNAseq_lib. Set RNASEQ_LIB_DIR.")
  }
}
for (f in c("plot_utils.R", "deg_utils.R", "enrichment_utils.R", "tcga_utils.R", "survival_utils.R")) {
  source(file.path(lib_dir, f))
}
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(survival); library(survminer); library(enrichplot)
})
theme_set(theme_publication())

# =============================================================================
# CONFIG — edit this block
# =============================================================================
RESULTS_RDATA <- "./1-DEG/TCGA_GEO_results.Rdata"
OUTDIR        <- "./3-Visualization/Custom"
DEG_PADJ_CUTOFF_OVERRIDE <- NULL   # NULL = reuse the value saved in the run
DEG_LFC_CUTOFF_OVERRIDE  <- NULL
TIME_UNIT_OVERRIDE       <- NULL   # NULL = reuse the run's TIME_UNIT

# Section switches
DO_VOLCANO     <- TRUE   # volcano for the Tumor_vs_Normal contrast
DO_BOXPLOT     <- TRUE   # single-gene expression boxplot by Tumor/Normal
DO_KM          <- TRUE   # KM by median + quartile for KMI genes
DO_COX_FOREST  <- TRUE   # forest plots from saved uni/multi Cox tables
DO_THEME       <- TRUE   # theme dot-heatmaps from saved ORA/GSEA objects

# Which genes to re-plot (boxplot + KM). NULL = reuse the run's GENES_FOR_SURVIVAL.
GENES_TO_PLOT  <- NULL

# =============================================================================
# Load saved results
# =============================================================================
if (!file.exists(RESULTS_RDATA)) stop("Not found: ", RESULTS_RDATA, " — run run_analysis.R first.")
load(RESULTS_RDATA)  # res_df, sig_genes, counts_raw, expr_log, clinical, surv_df,
                     # ego, ekegg, gsea_go, gsea_kegg, uni_cox, multi_cox,
                     # GENES_FOR_SURVIVAL, CLINICAL_VARS_FOR_KM, TIME_UNIT,
                     # SPECIES, TUMOR_NORMAL_DESIGN, DEG_PADJ_CUTOFF, DEG_LFC_CUTOFF

padj_cut  <- if (is.null(DEG_PADJ_CUTOFF_OVERRIDE)) DEG_PADJ_CUTOFF else DEG_PADJ_CUTOFF_OVERRIDE
lfc_cut   <- if (is.null(DEG_LFC_CUTOFF_OVERRIDE))  DEG_LFC_CUTOFF  else DEG_LFC_CUTOFF_OVERRIDE
time_unit <- if (is.null(TIME_UNIT_OVERRIDE))       TIME_UNIT       else TIME_UNIT_OVERRIDE
genes     <- if (is.null(GENES_TO_PLOT)) GENES_FOR_SURVIVAL else GENES_TO_PLOT
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
cat("Loaded results. Survival genes:", paste(genes, collapse = ", "), "\n")

# =============================================================================
# 1) Volcano plot
# =============================================================================
if (isTRUE(DO_VOLCANO)) {
  cat("\n[1] Volcano plot\n")
  if (is.null(res_df)) {
    message("  Skipped: no DEG result (TUMOR_NORMAL_DESIGN was FALSE).")
  } else {
    plot_volcano_pdf(res_df, comp_name = "Tumor_vs_Normal",
                     pvalue_thresh = padj_cut, log2fc_thresh = lfc_cut,
                     filename = file.path(OUTDIR, "Volcano_Tumor_vs_Normal.pdf"),
                     pvalue_column = "padj", lfc_column = "log2FoldChange")
    cat("  -> ", file.path(OUTDIR, "Volcano_Tumor_vs_Normal.pdf"), "\n")
  }
}

# =============================================================================
# 2) Single-gene expression boxplots
# =============================================================================
if (isTRUE(DO_BOXPLOT)) {
  cat("\n[2] Single-gene expression boxplots\n")
  if (!"condition" %in% colnames(clinical)) {
    message("  Skipped: clinical has no Tumor/Normal 'condition'.")
  } else {
    bdir <- file.path(OUTDIR, "Expression_boxplot"); dir.create(bdir, showWarnings = FALSE, recursive = TRUE)
    for (gene in genes) {
      if (!gene %in% rownames(expr_log)) { message("  Gene not in expression matrix: ", gene); next }
      plot_tcga_gene_boxplot_pdf(expr_log, gene, clinical$condition,
        filename = file.path(bdir, paste0("Expression_boxplot_", gene, ".pdf")),
        title = paste(gene, "expression"))
    }
    cat("  -> ", bdir, "\n")
  }
}

# =============================================================================
# 3) KM by median + quartile for chosen genes
# =============================================================================
if (isTRUE(DO_KM)) {
  cat("\n[3] KM curves (median + quartile)\n")
  if (is.null(surv_df)) {
    message("  Skipped: no survival data frame saved.")
  } else {
    kdir <- file.path(OUTDIR, "KM"); dir.create(kdir, showWarnings = FALSE, recursive = TRUE)
    for (gene in genes) {
      if (!gene %in% rownames(expr_log)) { message("  Gene not in expression matrix: ", gene); next }
      if (!gene %in% colnames(surv_df)) surv_df[[gene]] <- as.numeric(expr_log[gene, surv_df$barcode])
      plot_km_by_median_pdf(surv_df, value_col = gene,
        filename = file.path(kdir, paste0("KM_", gene, ".pdf")),
        title = paste(gene, "survival"), time_unit = time_unit)
      qcol <- paste0(gene, "_quartile")
      if (!qcol %in% colnames(surv_df)) surv_df[[qcol]] <- stratify_by_quantile(surv_df[[gene]], n_groups = 4)
      plot_km_by_group_pdf(surv_df, group_col = qcol,
        filename = file.path(kdir, paste0("KM_quartile_", gene, ".pdf")),
        title = paste(gene, "quartile survival"), time_unit = time_unit)
    }
    cat("  -> ", kdir, "\n")
  }
}

# =============================================================================
# 4) Cox forest plots (from saved Cox tables)
# =============================================================================
if (isTRUE(DO_COX_FOREST)) {
  cat("\n[4] Cox forest plots\n")
  cdir <- file.path(OUTDIR, "Cox_forest"); dir.create(cdir, showWarnings = FALSE, recursive = TRUE)
  if (!is.null(uni_cox) && nrow(uni_cox) > 0) {
    plot_cox_forest_pdf(uni_cox, filename = file.path(cdir, "univariate_Cox_forest.pdf"),
                        title = "Univariate Cox Regression")
  } else message("  No univariate Cox table saved.")
  if (!is.null(multi_cox) && nrow(multi_cox) > 0) {
    plot_cox_forest_pdf(multi_cox, filename = file.path(cdir, "multivariate_Cox_forest.pdf"),
                        title = "Multivariate Cox Regression")
  } else message("  No multivariate Cox table saved.")
  cat("  -> ", cdir, "\n")
}

# =============================================================================
# 5) Theme dot-heatmaps (from saved ORA/GSEA objects)
# =============================================================================
if (isTRUE(DO_THEME)) {
  cat("\n[5] Theme dot-heatmaps\n")
  theme_defs <- default_enrichment_themes()
  tdir <- file.path(OUTDIR, "ThemeEnrichment"); dir.create(tdir, showWarnings = FALSE, recursive = TRUE)
  go_ora_map    <- if (!is.null(ego))       list(Tumor_vs_Normal = ego)       else list()
  kegg_ora_map  <- if (!is.null(ekegg))     list(Tumor_vs_Normal = ekegg)     else list()
  go_gsea_map   <- if (!is.null(gsea_go))   list(Tumor_vs_Normal = gsea_go)   else list()
  kegg_gsea_map <- if (!is.null(gsea_kegg)) list(Tumor_vs_Normal = gsea_kegg) else list()
  if (length(go_ora_map) > 0)
    plot_theme_dotheatmap_from_results(go_ora_map, file.path(tdir, "Theme_dotheatmap_GO_ORA.pdf"),
      title = "GO ORA Biological Themes", subtitle = "GO-BP ORA", theme_defs = theme_defs, ontology_filter = "BP")
  if (length(kegg_ora_map) > 0)
    plot_theme_dotheatmap_from_results(kegg_ora_map, file.path(tdir, "Theme_dotheatmap_KEGG_ORA.pdf"),
      title = "KEGG ORA Pathway Themes", subtitle = "KEGG ORA", theme_defs = theme_defs, ontology_filter = NULL)
  if (length(go_gsea_map) > 0)
    plot_theme_dotheatmap_from_results(go_gsea_map, file.path(tdir, "Theme_dotheatmap_GO_GSEA.pdf"),
      title = "GO GSEA Biological Themes", subtitle = "GO-BP GSEA", theme_defs = theme_defs, ontology_filter = "BP")
  if (length(kegg_gsea_map) > 0)
    plot_theme_dotheatmap_from_results(kegg_gsea_map, file.path(tdir, "Theme_dotheatmap_KEGG_GSEA.pdf"),
      title = "KEGG GSEA Pathway Themes", subtitle = "KEGG GSEA", theme_defs = theme_defs, ontology_filter = NULL)
  if (length(c(go_ora_map, kegg_ora_map, go_gsea_map, kegg_gsea_map)) == 0)
    message("  No ORA/GSEA objects saved (TUMOR_NORMAL_DESIGN was FALSE).")
  cat("  -> ", tdir, "\n")
}

cat("\n========================================\n")
cat("Targeted visualization complete. Outputs under:", OUTDIR, "\n")
cat("========================================\n")
