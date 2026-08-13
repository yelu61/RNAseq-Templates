#!/usr/bin/env Rscript
# =============================================================================
# RNAseq_Limma_Voom — targeted re-visualization from saved results
# =============================================================================
# Run this AFTER run_analysis.R has produced results in the current project.
# It reads the saved intermediate objects — no voom / ORA / GSEA recompute:
#
#   1-DEG/limma_voom_results.Rdata   res_list, v, deg_summary, go_gsea_map, ...
#
# Use it to iterate cheaply on figures: restyle a volcano, pick specific GSEA
# terms, or rebuild a theme dot-heatmap — without re-running the pipeline.
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
for (f in c("plot_utils.R", "io_utils.R", "data_utils.R", "deg_utils.R", "enrichment_utils.R")) {
  source(file.path(lib_dir, f))
}
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(enrichplot)
})
theme_set(theme_publication())

# =============================================================================
# CONFIG — edit this block
# =============================================================================
RESULTS_RDATA <- "./1-DEG/limma_voom_results.Rdata"
OUTDIR        <- "./3-Visualization/Custom"
DEG_PADJ_CUTOFF_OVERRIDE <- NULL   # NULL = reuse the value saved in the run
DEG_LFC_CUTOFF_OVERRIDE  <- NULL

# Section switches
DO_VOLCANO     <- TRUE   # volcano per comparison
DO_HEATMAP     <- TRUE   # top-DEG expression heatmap (needs >= 2 sig genes)
DO_GSEA_TERMS  <- TRUE   # one gseaplot2 figure per chosen term
DO_THEME       <- TRUE   # theme dot-heatmaps from saved ORA/GSEA maps

# GSEA single-term selection
GSEA_COMPARISON <- NULL                 # NULL = first comparison
GSEA_DB         <- "GO"                 # "GO" or "KEGG"
GSEA_TERM_MATCH <- "top"                # character vector of Description substrings, or "top"
GSEA_TOP_N      <- 6                    # used when GSEA_TERM_MATCH = "top"

# =============================================================================
# Load saved results
# =============================================================================
if (!file.exists(RESULTS_RDATA)) stop("Not found: ", RESULTS_RDATA, " — run run_analysis.R first.")
load(RESULTS_RDATA)  # res_list, deg_summary, v, countData, GROUP_LEVELS, COMPARISONS, SPECIES,
                     # DEG_PADJ_CUTOFF, DEG_LFC_CUTOFF, go_ora_map, kegg_ora_map, go_gsea_map, kegg_gsea_map

padj_cut <- if (is.null(DEG_PADJ_CUTOFF_OVERRIDE)) DEG_PADJ_CUTOFF else DEG_PADJ_CUTOFF_OVERRIDE
lfc_cut  <- if (is.null(DEG_LFC_CUTOFF_OVERRIDE))  DEG_LFC_CUTOFF  else DEG_LFC_CUTOFF_OVERRIDE
group_colors <- make_group_colors(GROUP_LEVELS)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
cat("Loaded results. Comparisons:", paste(names(res_list), collapse = ", "), "\n")

# =============================================================================
# 1) Volcano plots
# =============================================================================
if (isTRUE(DO_VOLCANO)) {
  cat("\n[1] Volcano plots\n")
  vdir <- file.path(OUTDIR, "Volcano"); dir.create(vdir, showWarnings = FALSE, recursive = TRUE)
  for (comp_name in names(res_list)) {
    plot_volcano_pdf(res_list[[comp_name]], comp_name = comp_name,
                     pvalue_thresh = padj_cut, log2fc_thresh = lfc_cut,
                     filename = file.path(vdir, paste0("Volcano_", comp_name, ".pdf")),
                     pvalue_column = "padj", lfc_column = "log2FoldChange")
  }
  cat("  -> ", vdir, "\n")
}

# =============================================================================
# 2) Top-DEG heatmap
# =============================================================================
if (isTRUE(DO_HEATMAP)) {
  cat("\n[2] Top-DEG heatmap\n")
  sig <- unique(unlist(lapply(res_list, function(res) {
    res |>
      dplyr::filter(!is.na(padj), padj < padj_cut, abs(log2FoldChange) > lfc_cut) |>
      dplyr::pull(gene_name)
  })))
  if (length(sig) < 2) {
    message("  Skipped: fewer than 2 significant genes at the chosen thresholds.")
  } else {
    mat <- v$E[intersect(sig, rownames(v$E)), ]
    plot_expression_heatmap_pdf(mat,
      filename = file.path(OUTDIR, "DEG_heatmap.pdf"), title = "limma-voom DEGs",
      group = GROUPS, group_levels = GROUP_LEVELS, group_colors = group_colors,
      width = mm_to_in(183), height = mm_to_in(247))
    cat("  -> ", file.path(OUTDIR, "DEG_heatmap.pdf"), "\n")
  }
}

# =============================================================================
# 3) GSEA single-term figures (reads cached gseaResult objects)
# =============================================================================
if (isTRUE(DO_GSEA_TERMS)) {
  cat("\n[3] GSEA single-term figures\n")
  gsea_map <- if (toupper(GSEA_DB) == "KEGG") kegg_gsea_map else go_gsea_map
  comp <- if (!is.null(GSEA_COMPARISON)) GSEA_COMPARISON else names(gsea_map)[1]
  gsea_obj <- gsea_map[[comp]]
  if (is.null(gsea_obj) || nrow(as.data.frame(gsea_obj)) == 0) {
    message("  Skipped: no ", GSEA_DB, " gseaResult for comparison ", comp)
  } else {
    gdf <- as.data.frame(gsea_obj); gdf <- gdf[!is.na(gdf$NES), , drop = FALSE]
    if (length(GSEA_TERM_MATCH) == 1 && identical(tolower(GSEA_TERM_MATCH), "top")) {
      gdf <- gdf[order(gdf$p.adjust, -abs(gdf$NES)), , drop = FALSE]
      picked <- utils::head(gdf, GSEA_TOP_N)
    } else {
      pat <- paste(GSEA_TERM_MATCH, collapse = "|")
      picked <- gdf[grepl(pat, gdf$Description, ignore.case = TRUE), , drop = FALSE]
    }
    if (nrow(picked) == 0) {
      message("  No GSEA terms matched.")
    } else {
      term_dir <- file.path(OUTDIR, "GSEA_terms", paste0(GSEA_DB, "_", comp))
      dir.create(term_dir, showWarnings = FALSE, recursive = TRUE)
      plot_gsea_term_figures_from_df(gsea_obj, picked, outdir = term_dir,
                                     contrast_label = comp, prefix = paste0("gseaplot2_", GSEA_DB))
      cat("  ", nrow(picked), "term figure(s) ->", term_dir, "\n")
    }
  }
}

# =============================================================================
# 4) Theme dot-heatmaps (from saved ORA/GSEA maps)
# =============================================================================
if (isTRUE(DO_THEME)) {
  cat("\n[4] Theme dot-heatmaps\n")
  theme_defs <- default_enrichment_themes()
  tdir <- file.path(OUTDIR, "ThemeEnrichment"); dir.create(tdir, showWarnings = FALSE, recursive = TRUE)
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
  cat("  -> ", tdir, "\n")
}

cat("\n========================================\n")
cat("Targeted visualization complete. Outputs under:", OUTDIR, "\n")
cat("========================================\n")
