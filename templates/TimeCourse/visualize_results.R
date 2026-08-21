#!/usr/bin/env Rscript
# =============================================================================
# RNAseq_TimeCourse — targeted re-visualization from saved results
# =============================================================================
# Run this AFTER run_analysis.R has produced results in the current project.
# It reads the saved intermediate objects — no Mfuzz / DESeq2 recompute:
#
#   5-TimeCourse/timecourse_results.Rdata   expr_mean, cluster_df, mfuzz_result,
#                                           eset, ora_results, tp_res_list, ...
#
# Use it to iterate cheaply on figures: restyle an Mfuzz trend/heatmap, rebuild a
# theme dot-heatmap, or re-cut a time-point DEG volcano — without re-running the
# pipeline.
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
for (f in c("plot_utils.R", "io_utils.R", "data_utils.R", "deg_utils.R",
            "enrichment_utils.R", "timecourse_utils.R")) {
  source(file.path(lib_dir, f))
}
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})
theme_set(theme_publication())

# =============================================================================
# CONFIG — edit this block
# =============================================================================
RESULTS_RDATA <- "./5-TimeCourse/timecourse_results.Rdata"
OUTDIR        <- "./3-Visualization/Custom"
DEG_PADJ_CUTOFF_OVERRIDE <- NULL   # NULL = reuse the value saved in the run
DEG_LFC_CUTOFF_OVERRIDE  <- NULL

# Section switches
DO_MFUZZ_TRENDS  <- TRUE   # Mfuzz cluster trend panels (needs saved eset + mfuzz_result)
DO_MFUZZ_HEATMAP <- TRUE   # core-gene time-course heatmap (needs >= 2 core genes)
DO_THEME         <- TRUE   # theme dot-heatmap from the saved per-cluster ORA map
DO_TIMEPLOTS     <- TRUE   # time-point DEG volcano + summary from saved tp_res_list

# Mfuzz trend restyle
TREND_WIDTH  <- 14
TREND_HEIGHT <- 10

# =============================================================================
# Load saved results
# =============================================================================
if (!file.exists(RESULTS_RDATA)) stop("Not found: ", RESULTS_RDATA, " — run run_analysis.R first.")
load(RESULTS_RDATA)  # expr_mean, cluster_df, mfuzz_result, eset, ora_results,
                     # tp_res_list, tp_summary, TIME_LEVELS, SPECIES, RUN_MFUZZ,
                     # RUN_TIMEPOINT_DEG, baseline_time, DEG_PADJ_CUTOFF, DEG_LFC_CUTOFF

padj_cut <- if (is.null(DEG_PADJ_CUTOFF_OVERRIDE)) DEG_PADJ_CUTOFF else DEG_PADJ_CUTOFF_OVERRIDE
lfc_cut  <- if (is.null(DEG_LFC_CUTOFF_OVERRIDE))  DEG_LFC_CUTOFF  else DEG_LFC_CUTOFF_OVERRIDE
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
cat("Loaded results. Time points:", paste(TIME_LEVELS, collapse = ", "),
    "| Mfuzz:", isTRUE(RUN_MFUZZ), "| timepoint DEG:", isTRUE(RUN_TIMEPOINT_DEG), "\n")

# =============================================================================
# 1) Mfuzz cluster trend panels
# =============================================================================
if (isTRUE(DO_MFUZZ_TRENDS)) {
  cat("\n[1] Mfuzz cluster trends\n")
  if (!isTRUE(RUN_MFUZZ) || is.null(eset) || is.null(mfuzz_result)) {
    message("  Skipped: no saved Mfuzz objects (RUN_MFUZZ was FALSE or clustering was skipped).")
  } else if (!requireNamespace("Mfuzz", quietly = TRUE)) {
    message("  Skipped: package 'Mfuzz' is not installed in this session.")
  } else {
    # mfuzz.plot() calls Biobase::exprs() unqualified, so attach Mfuzz (which
    # brings Biobase onto the search path) rather than relying on requireNamespace.
    suppressPackageStartupMessages(library(Mfuzz))
    plot_mfuzz_trends_pdf(
      eset, mfuzz_result,
      filename = file.path(OUTDIR, "mfuzz_trends.pdf"),
      time_labels = TIME_LEVELS,
      width = TREND_WIDTH, height = TREND_HEIGHT
    )
    cat("  -> ", file.path(OUTDIR, "mfuzz_trends.pdf"), "\n")
  }
}

# =============================================================================
# 2) Mfuzz core-gene heatmap (re-plotted from the mean-by-time matrix)
# =============================================================================
if (isTRUE(DO_MFUZZ_HEATMAP)) {
  cat("\n[2] Mfuzz core-gene heatmap\n")
  core_df <- if (!is.null(cluster_df)) cluster_df[cluster_df$core_gene, ] else NULL
  if (!isTRUE(RUN_MFUZZ) || is.null(core_df)) {
    message("  Skipped: no saved cluster table.")
  } else if (nrow(core_df) < 2) {
    message("  Skipped: fewer than 2 core genes.")
  } else if (!(requireNamespace("ComplexHeatmap", quietly = TRUE) && requireNamespace("circlize", quietly = TRUE))) {
    message("  Skipped: ComplexHeatmap/circlize not installed.")
  } else {
    group_colors <- make_group_colors(TIME_LEVELS)
    # Re-plot against the per-time-point mean matrix (one column per time point).
    plot_timecourse_heatmap_pdf(
      expr_mean, core_df,
      group_vec = TIME_LEVELS,
      group_levels = TIME_LEVELS,
      group_colors = group_colors,
      filename = file.path(OUTDIR, "mfuzz_core_heatmap.pdf"),
      width = mm_to_in(183), height = mm_to_in(247)
    )
    cat("  -> ", file.path(OUTDIR, "mfuzz_core_heatmap.pdf"), "\n")
  }
}

# =============================================================================
# 3) Theme dot-heatmap (from the saved per-cluster ORA map)
# =============================================================================
if (isTRUE(DO_THEME)) {
  cat("\n[3] Theme dot-heatmap\n")
  # vapply (not the notebook's sapply) so an empty ora_results does not crash.
  drop_empty <- function(lst) lst[vapply(lst, function(x) !is.null(x) && nrow(as.data.frame(x)) > 0, logical(1))]
  ora_map <- drop_empty(ora_results)
  if (length(ora_map) == 0) {
    message("  Skipped: no non-empty per-cluster ORA results saved.")
  } else {
    theme_defs <- default_enrichment_themes()
    tdir <- file.path(OUTDIR, "ThemeEnrichment"); dir.create(tdir, showWarnings = FALSE, recursive = TRUE)
    plot_theme_dotheatmap_from_results(
      ora_map,
      filename = file.path(tdir, "Theme_dotheatmap_GO_ORA_mfuzz_clusters.pdf"),
      title = "GO ORA Biological Themes (Mfuzz Clusters)",
      subtitle = "GO-BP ORA | top terms per theme per cluster",
      theme_defs = theme_defs, ontology_filter = "BP", top_n = 6
    )
    cat("  -> ", tdir, "\n")
  }
}

# =============================================================================
# 4) Time-point DEG volcano + summary (from saved tp_res_list)
# =============================================================================
if (isTRUE(DO_TIMEPLOTS)) {
  cat("\n[4] Time-point DEG volcano + summary\n")
  if (!isTRUE(RUN_TIMEPOINT_DEG) || length(tp_res_list) == 0) {
    message("  Skipped: no saved time-point DEG results.")
  } else {
    vdir <- file.path(OUTDIR, "Volcano"); dir.create(vdir, showWarnings = FALSE, recursive = TRUE)
    for (comp_name in names(tp_res_list)) {
      plot_volcano_pdf(
        tp_res_list[[comp_name]], comp_name = comp_name,
        pvalue_thresh = padj_cut, log2fc_thresh = lfc_cut,
        pvalue_column = "padj", lfc_column = "log2FoldChange_shrunken",
        filename = file.path(vdir, paste0("Volcano_", comp_name, ".pdf"))
      )
    }
    cat("  -> ", vdir, "\n")
    if (!is.null(tp_summary) && nrow(tp_summary) > 0) {
      plot_timepoint_deg_summary_pdf(
        tp_summary,
        filename = file.path(OUTDIR, "Timepoint_DEG_summary.pdf")
      )
      cat("  -> ", file.path(OUTDIR, "Timepoint_DEG_summary.pdf"), "\n")
    }
  }
}

cat("\n========================================\n")
cat("Targeted visualization complete. Outputs under:", OUTDIR, "\n")
cat("========================================\n")
