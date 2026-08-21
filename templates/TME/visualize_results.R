#!/usr/bin/env Rscript
# =============================================================================
# RNAseq_TME — targeted re-visualization from saved results
# =============================================================================
# Run this AFTER run_analysis.R has produced results in the current project.
# It reads the saved intermediate objects — no deconvolution / ssGSEA recompute:
#
#   4-TME/tme_results.Rdata   iobr_results, estimate_scores, native_cibersort,
#                             ssgsea_scores, meta, group_colors, ...
#
# Use it to iterate cheaply on figures: restyle a barplot/boxplot, rebuild an
# ESTIMATE boxplot or an ssGSEA heatmap — without re-running the pipeline.
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
file_arg <- gsub("~\\+~", " ", file_arg)
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
for (f in c("plot_utils.R", "tme_utils.R", "data_utils.R")) {
  source(file.path(lib_dir, f))
}
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(pheatmap)
})
theme_set(theme_publication())

# =============================================================================
# CONFIG — edit this block
# =============================================================================
RESULTS_RDATA <- "./4-TME/tme_results.Rdata"
OUTDIR        <- "./4-TME/Custom"

# Section switches
DO_TME_BARPLOT   <- TRUE   # stacked barplot + boxplot for a fraction method
DO_TME_HEATMAP   <- TRUE   # heatmap for a score/fraction method
DO_ESTIMATE      <- TRUE   # ESTIMATE scores boxplot (IOBR estimate preferred)
DO_SSGSEA_HEATMAP <- TRUE  # ssGSEA signature heatmap

# Method selection for the barplot/boxplot and heatmap sections. These refer to
# names(iobr_results), e.g. "cibersort", "epic", "xcell", "estimate".
BARPLOT_METHOD  <- "cibersort"   # fraction-based method for barplot + boxplot
HEATMAP_METHOD  <- "xcell"       # method for the TME heatmap

# =============================================================================
# Load saved results
# =============================================================================
if (!file.exists(RESULTS_RDATA)) stop("Not found: ", RESULTS_RDATA, " — run run_analysis.R first.")
load(RESULTS_RDATA)  # expr_tme, iobr_results, tme_combined, estimate_scores, native_cibersort,
                     # ssgsea_scores, gs, meta, group_df_for_plot, group_colors,
                     # SAMPLE_COLUMN, GROUP_COLUMN, GROUP_LEVELS, SPECIES, RUN_* flags

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
cat("Loaded results. IOBR methods:", paste(names(iobr_results), collapse = ", "), "\n")

# =============================================================================
# 1) TME stacked barplot + boxplot (fraction-based method)
# =============================================================================
if (isTRUE(DO_TME_BARPLOT)) {
  cat("\n[1] TME barplot + boxplot (", BARPLOT_METHOD, ")\n")
  if (!BARPLOT_METHOD %in% names(iobr_results)) {
    message("  Skipped: method '", BARPLOT_METHOD, "' not present in iobr_results.")
  } else {
    method_label <- toupper(BARPLOT_METHOD)
    method_long <- melt_tme_results(iobr_results[[BARPLOT_METHOD]], id_column = SAMPLE_COLUMN,
      group_df = group_df_for_plot, sample_col = SAMPLE_COLUMN, group_col = GROUP_COLUMN) |>
      dplyr::filter(!grepl("P-value|Correlation|RMSE", .data$cell_type))
    bar_size <- calc_tme_barplot_size(length(unique(method_long[[SAMPLE_COLUMN]])), length(unique(method_long$cell_type)))
    box_size <- calc_tme_boxplot_size(length(unique(method_long$cell_type)))
    plot_tme_barplot_pdf(method_long, group_col = GROUP_COLUMN, sample_col = SAMPLE_COLUMN,
      filename = file.path(OUTDIR, paste0("IOBR_", method_label, "_barplot.pdf")),
      title = paste(method_label, "Cell Fractions"),
      width = bar_size["width"], height = bar_size["height"])
    plot_tme_boxplot_pdf(method_long, group_col = GROUP_COLUMN, value_col = "fraction",
      filename = file.path(OUTDIR, paste0("IOBR_", method_label, "_boxplot.pdf")),
      title = paste(method_label, "Cell Fractions by Group"),
      width = box_size["width"], height = box_size["height"], group_colors = group_colors)
    cat("  -> ", OUTDIR, "\n")
  }
}

# =============================================================================
# 2) TME heatmap (score/fraction method)
# =============================================================================
if (isTRUE(DO_TME_HEATMAP)) {
  cat("\n[2] TME heatmap (", HEATMAP_METHOD, ")\n")
  if (!HEATMAP_METHOD %in% names(iobr_results)) {
    message("  Skipped: method '", HEATMAP_METHOD, "' not present in iobr_results.")
  } else {
    plot_tme_heatmap_pdf(iobr_results[[HEATMAP_METHOD]], meta, group_col = GROUP_COLUMN,
      sample_col = SAMPLE_COLUMN, group_colors = group_colors,
      filename = file.path(OUTDIR, paste0("IOBR_", toupper(HEATMAP_METHOD), "_heatmap.pdf")),
      title = paste(toupper(HEATMAP_METHOD), "Scores"), width = 10, height = 12)
    cat("  -> ", file.path(OUTDIR, paste0("IOBR_", toupper(HEATMAP_METHOD), "_heatmap.pdf")), "\n")
  }
}

# =============================================================================
# 3) ESTIMATE scores boxplot
# =============================================================================
if (isTRUE(DO_ESTIMATE)) {
  cat("\n[3] ESTIMATE scores boxplot\n")
  est_source <- if ("estimate" %in% names(iobr_results)) iobr_results[["estimate"]] else estimate_scores
  if (is.null(est_source)) {
    message("  Skipped: no ESTIMATE scores in the saved results (neither IOBR estimate nor native).")
  } else {
    est_long <- melt_estimate_scores(est_source, id_column = SAMPLE_COLUMN,
                                     group_df = group_df_for_plot,
                                     sample_col = SAMPLE_COLUMN, group_col = GROUP_COLUMN)
    if (is.null(est_long) || nrow(est_long) == 0) {
      message("  Skipped: ESTIMATE table has no recognized score columns.")
    } else {
      plot_estimate_boxplot_pdf(est_long, group_col = GROUP_COLUMN,
        filename = file.path(OUTDIR, "ESTIMATE_scores_boxplot.pdf"),
        title = "ESTIMATE Scores by Group", group_colors = group_colors)
      cat("  -> ", file.path(OUTDIR, "ESTIMATE_scores_boxplot.pdf"), "\n")
    }
  }
}

# =============================================================================
# 4) ssGSEA signature heatmap
# =============================================================================
if (isTRUE(DO_SSGSEA_HEATMAP)) {
  cat("\n[4] ssGSEA signature heatmap\n")
  if (is.null(ssgsea_scores)) {
    message("  Skipped: no ssGSEA scores in the saved results.")
  } else {
    ann <- data.frame(Group = as.character(meta[[GROUP_COLUMN]]))
    rownames(ann) <- meta[[SAMPLE_COLUMN]]
    ann$Group <- factor(ann$Group, levels = GROUP_LEVELS)
    annotation_colors <- list(Group = group_colors)
    pheatmap(ssgsea_scores, annotation_col = ann, annotation_colors = annotation_colors,
             scale = "row", cluster_cols = TRUE, cluster_rows = TRUE,
             filename = file.path(OUTDIR, "ssGSEA_heatmap.pdf"), width = 8, height = 7)
    cat("  -> ", file.path(OUTDIR, "ssGSEA_heatmap.pdf"), "\n")
  }
}

cat("\n========================================\n")
cat("Targeted visualization complete. Outputs under:", OUTDIR, "\n")
cat("========================================\n")
