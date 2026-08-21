#!/usr/bin/env Rscript
# =============================================================================
# RNAseq_General — targeted re-visualization from saved results
# =============================================================================
# Run this AFTER run_analysis.R has produced results in the current project.
# It reads the saved intermediate objects — no DESeq2 / ORA / GSEA recompute:
#
#   1-DEG/DEG_results.Rdata   res_list, vsd_mat, colData, GROUP_LEVELS, ...
#   2-GSEA/gsea_results.rds   cached gseaResult objects (for gseaplot2 figures)
#   2-GSEA/<threshold>/*.csv  ORA result tables (for theme dot-heatmaps)
#
# Use it to iterate cheaply on figures: change key genes, pick specific GSEA
# terms, or restyle a theme dot-heatmap — without re-running the pipeline.
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
for (f in c("plot_utils.R", "io_utils.R", "data_utils.R", "deg_utils.R", "enrichment_utils.R")) {
  source(file.path(lib_dir, f))
}
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(stringr)
  library(enrichplot); library(DESeq2)
})
theme_set(theme_publication())

# =============================================================================
# CONFIG — edit this block
# =============================================================================
DEG_RDATA   <- "./1-DEG/DEG_results.Rdata"
GSEA_RDS    <- "./2-GSEA/gsea_results.rds"
GSEA_DIR    <- "./2-GSEA"                 # where GSEA_*_<comp>.csv live
ORA_DIR     <- "./2-GSEA"                 # ORA csvs live under <threshold>/
OUTDIR      <- "./3-Visualization/Custom"
DEG_PVALUE_COLUMN <- "padj"
DEG_LFC_COLUMN    <- "log2FoldChange_raw"
DEFAULT_THRESHOLD <- "standard"

# Section switches
DO_KEY_GENES      <- TRUE    # single-gene bar/SEM + key-gene heatmap
DO_GSEA_TERMS     <- TRUE    # one gseaplot2 figure per chosen term (needs gsea_results.rds)
DO_ORA_THEME      <- TRUE    # ORA theme dot-heatmap from saved ORA csvs

# 1) Key genes ---------------------------------------------------------------
KEY_GENES <- c("Tnf", "Il1b", "Il6", "Cxcl10", "Nos2")

# 2) GSEA single-term figures -------------------------------------------------
# Which comparison(s) and which terms to draw. Use a character vector of term
# Description substrings (matched case-insensitively against that comparison's
# GSEA table), or "top" to take the top N by |NES|.
GSEA_COMPARISON <- NULL                    # NULL = first comparison found
GSEA_TERM_MATCH <- c("inflammatory", "chemokine", "cytokine")  # substrings, OR
GSEA_TOP_N      <- 6                       # used when GSEA_TERM_MATCH = "top"
GSEA_DB         <- "GO"                    # "GO" or "KEGG"

# 3) ORA theme dot-heatmap ----------------------------------------------------
ORA_THEME_DB      <- "GO"                  # "GO" or "KEGG"
ORA_THEME_ONTOLOGY <- "BP"                 # GO only; NULL for KEGG
ORA_THEME_TOP_N   <- 6
ORA_THEME_DEFS    <- default_enrichment_themes()  # supply your own to customize

# =============================================================================
# Load saved results
# =============================================================================
if (!file.exists(DEG_RDATA)) stop("Not found: ", DEG_RDATA, " — run run_analysis.R first.")
load(DEG_RDATA)  # res_list, deg_by_threshold, dds, vsd, vsd_mat, colData, GROUP_LEVELS, COMPARISONS, SPECIES, ...

if (!exists("GROUP_LEVELS") || is.null(GROUP_LEVELS)) {
  GROUP_LEVELS <- levels(factor(colData$condition))
}
group_colors <- make_group_colors(GROUP_LEVELS)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
cat("Loaded results. Comparisons:", paste(names(res_list), collapse = ", "), "\n")

# =============================================================================
# 1) Key genes: single-gene bar/SEM + heatmap
# =============================================================================
if (isTRUE(DO_KEY_GENES)) {
  cat("\n[1] Key-gene visualization\n")
  kg_dir <- file.path(OUTDIR, "KeyGenes")
  dir.create(kg_dir, showWarnings = FALSE, recursive = TRUE)

  plot_gene <- function(gene) {
    if (!gene %in% rownames(vsd_mat)) { cat("  not found:", gene, "\n"); return(invisible(NULL)) }
    pd <- data.frame(
      sample = colnames(vsd_mat),
      expression = as.numeric(vsd_mat[gene, ]),
      condition = factor(as.character(colData[colnames(vsd_mat), "condition"]), levels = GROUP_LEVELS))
    plot_group_bar_sem_pdf(pd, value_col = "expression", group_col = "condition",
                           comparisons = combn(GROUP_LEVELS, 2, simplify = FALSE),
                           method = "t.test", p_adjust_method = "BH",
                           title = gene, ylab = "VST Expression", group_colors = group_colors,
                           filename = file.path(kg_dir, paste0("Barplot_", gene, ".pdf")),
                           width = 5.5, height = 6, y_from_zero = TRUE)
  }
  invisible(lapply(KEY_GENES, plot_gene))

  found <- intersect(KEY_GENES, rownames(vsd_mat))
  if (length(found) > 0) {
    plot_expression_heatmap_pdf(vsd_mat[found, , drop = FALSE],
                                filename = file.path(kg_dir, "Heatmap_key_genes.pdf"),
                                title = "Key Genes",
                                group = as.character(colData[colnames(vsd_mat), "condition"]),
                                group_levels = GROUP_LEVELS, group_colors = group_colors,
                                show_row_names = TRUE, show_column_names = TRUE,
                                cluster_rows = TRUE, cluster_columns = TRUE, row_font_size = 8,
                                width = 8, height = max(5, 0.28 * length(found) + 3))
  }
  cat("  Key-gene figures ->", kg_dir, "\n")
}

# =============================================================================
# 2) GSEA single-term figures (reads cached gseaResult objects)
# =============================================================================
if (isTRUE(DO_GSEA_TERMS)) {
  cat("\n[2] GSEA single-term figures\n")
  if (!file.exists(GSEA_RDS)) {
    message("  Skipped: ", GSEA_RDS, " not found. Re-run run_analysis.R (it now caches gseaResult objects),")
    message("  or set DO_GSEA_TERMS <- FALSE.")
  } else {
    gsea_results <- readRDS(GSEA_RDS)
    comp <- if (!is.null(GSEA_COMPARISON)) GSEA_COMPARISON else names(gsea_results)[1]
    if (!comp %in% names(gsea_results)) {
      message("  Comparison '", comp, "' not in gsea_results; available: ", paste(names(gsea_results), collapse = ", "))
    } else {
      gsea_obj <- if (toupper(GSEA_DB) == "KEGG") gsea_results[[comp]]$kegg else gsea_results[[comp]]$go
      if (is.null(gsea_obj) || nrow(as.data.frame(gsea_obj)) == 0) {
        message("  No ", GSEA_DB, " gseaResult for comparison ", comp)
      } else {
        gdf <- as.data.frame(gsea_obj)
        gdf <- gdf[!is.na(gdf$NES), , drop = FALSE]
        # pick terms
        if (length(GSEA_TERM_MATCH) == 1 && identical(tolower(GSEA_TERM_MATCH), "top")) {
          gdf <- gdf[order(gdf$p.adjust, -abs(gdf$NES)), , drop = FALSE]
          picked <- utils::head(gdf, GSEA_TOP_N)
        } else {
          pat <- paste(GSEA_TERM_MATCH, collapse = "|")
          picked <- gdf[grepl(pat, gdf$Description, ignore.case = TRUE), , drop = FALSE]
          picked <- picked[order(picked$p.adjust, -abs(picked$NES)), , drop = FALSE]
        }
        if (nrow(picked) == 0) {
          message("  No GSEA terms matched. Check GSEA_TERM_MATCH against the Description column.")
        } else {
          term_dir <- file.path(OUTDIR, "GSEA_terms", paste0(GSEA_DB, "_", comp))
          dir.create(term_dir, showWarnings = FALSE, recursive = TRUE)
          plot_gsea_term_figures_from_df(gsea_obj, picked, outdir = term_dir,
                                         contrast_label = comp, prefix = paste0("gseaplot2_", GSEA_DB))
          cat("  ", nrow(picked), "term figure(s) ->", term_dir, "\n")
        }
      }
    }
  }
}

# =============================================================================
# 3) ORA theme dot-heatmap (reads saved ORA csvs)
# =============================================================================
if (isTRUE(DO_ORA_THEME)) {
  cat("\n[3] ORA theme dot-heatmap\n")
  th_dir <- file.path(ORA_DIR, DEFAULT_THRESHOLD)
  prefix <- paste0(ORA_THEME_DB, "_ORA_")
  ora_files <- list.files(th_dir, pattern = paste0("^", prefix, ".*\\.csv$"), full.names = TRUE)
  # exclude the UP/DOWN-split tables so each comparison appears once
  ora_files <- ora_files[!grepl("_UP_|_DOWN_|_bidirectional_", basename(ora_files))]
  if (length(ora_files) == 0) {
    message("  Skipped: no ", prefix, "*.csv under ", th_dir)
  } else {
    result_map <- list()
    for (f in ora_files) {
      comp <- sub(paste0("^", prefix), "", tools::file_path_sans_ext(basename(f)))
      df <- utils::read.csv(f, stringsAsFactors = FALSE)
      if (nrow(df) > 0) result_map[[comp]] <- df
    }
    if (length(result_map) == 0) {
      message("  Skipped: ORA csvs were empty.")
    } else {
      ont <- if (toupper(ORA_THEME_DB) == "KEGG") NULL else ORA_THEME_ONTOLOGY
      theme_out <- file.path(OUTDIR, "ORA_theme")
      dir.create(theme_out, showWarnings = FALSE, recursive = TRUE)
      p <- plot_theme_dotheatmap_from_results(
        result_map,
        filename = file.path(theme_out, paste0("Theme_dotheatmap_", ORA_THEME_DB, "_ORA.pdf")),
        title = paste(ORA_THEME_DB, "ORA Biological Themes"),
        subtitle = paste0(ORA_THEME_DB, if (!is.null(ont)) paste0("-", ont) else "", " ORA | ", DEFAULT_THRESHOLD),
        theme_defs = ORA_THEME_DEFS, ontology_filter = ont, top_n = ORA_THEME_TOP_N)
      if (!is.null(p)) cat("  Theme dot-heatmap ->", theme_out, "\n")
    }
  }
}

cat("\n========================================\n")
cat("Targeted visualization complete. Outputs under:", OUTDIR, "\n")
cat("========================================\n")
