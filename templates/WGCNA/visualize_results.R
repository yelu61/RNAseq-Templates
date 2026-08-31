#!/usr/bin/env Rscript
# =============================================================================
# RNAseq_WGCNA — targeted re-visualization from saved results
# =============================================================================
# Run this AFTER run_analysis.R has produced results in the current project.
# It reads the saved intermediate objects — no blockwiseModules recompute:
#
#   5-WGCNA/WGCNA_results.Rdata   net, MEs, moduleColors, datExpr, traits,
#                                 moduleTraitCor, moduleTraitP, sft, ...
#
# Use it to iterate cheaply on figures: restyle the module-trait heatmap,
# re-export hub-gene tables for chosen modules, or rebuild the dendrogram and
# soft-threshold plots — without re-running the network construction.
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
for (f in c("plot_utils.R", "data_utils.R")) {
  source(file.path(lib_dir, f))
}
suppressPackageStartupMessages({
  library(WGCNA); library(ggplot2); library(dplyr)
})
theme_set(theme_publication())

# =============================================================================
# CONFIG — edit this block
# =============================================================================
RESULTS_RDATA <- "./5-WGCNA/WGCNA_results.Rdata"
OUTDIR        <- "./5-WGCNA/Custom"

# Section switches
DO_SOFT_THRESHOLD <- TRUE   # rebuild the soft-threshold diagnostic plot
DO_DENDROGRAM     <- TRUE   # module dendrogram + colors
DO_MT_HEATMAP     <- TRUE   # restyle the module-trait heatmap
DO_HUB            <- TRUE   # per-module hub-gene tables

# Hub-gene export controls
HUB_MODULES <- NULL         # NULL = all non-grey modules, or e.g. c("blue", "turquoise")
HUB_TOP_N   <- NULL         # NULL = full table, or keep the top-N genes by |kME|

# =============================================================================
# Load saved results
# =============================================================================
if (!file.exists(RESULTS_RDATA)) stop("Not found: ", RESULTS_RDATA, " — run run_analysis.R first.")
load(RESULTS_RDATA)  # net, MEs, moduleColors, datExpr, traits, trait_numeric, moduleTraitCor,
                     # moduleTraitP, sft, soft_power, all_modules, hub_list, NETWORK_TYPE, ...

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
cat("Loaded results. Modules:", paste(all_modules, collapse = ", "),
    "| soft power =", soft_power, "\n")

# =============================================================================
# 1) Soft-threshold diagnostic plot
# =============================================================================
if (isTRUE(DO_SOFT_THRESHOLD)) {
  cat("\n[1] Soft-threshold selection plot\n")
  if (is.null(sft) || is.null(sft$fitIndices)) {
    message("  Skipped: no 'sft' object in the saved results.")
  } else {
    plot_wgcna_soft_threshold_pdf(
      sft, selected_power = soft_power,
      filename = file.path(OUTDIR, "Soft_threshold_selection.pdf")
    )
    cat("  -> ", file.path(OUTDIR, "Soft_threshold_selection.pdf"), "\n")
  }
}

# =============================================================================
# 2) Module dendrogram
# =============================================================================
if (isTRUE(DO_DENDROGRAM)) {
  cat("\n[2] Module dendrogram\n")
  if (is.null(net) || is.null(moduleColors)) {
    message("  Skipped: no 'net'/'moduleColors' in the saved results.")
  } else {
    plot_wgcna_module_dendrogram_pdf(
      net, moduleColors, file.path(OUTDIR, "Module_dendrogram.pdf")
    )
    cat("  -> ", file.path(OUTDIR, "Module_dendrogram.pdf"), "\n")
  }
}

# =============================================================================
# 3) Module-trait heatmap (restyle)
# =============================================================================
if (isTRUE(DO_MT_HEATMAP)) {
  cat("\n[3] Module-trait heatmap\n")
  if (is.null(moduleTraitCor) || is.null(moduleTraitP) || is.null(trait_numeric)) {
    message("  Skipped: module-trait matrices not found in the saved results.")
  } else {
    plot_wgcna_module_trait_heatmap_pdf(
      moduleTraitCor, moduleTraitP, file.path(OUTDIR, "Module_trait_heatmap.pdf")
    )
    cat("  -> ", file.path(OUTDIR, "Module_trait_heatmap.pdf"), "\n")
  }
}

# =============================================================================
# 4) Per-module hub-gene export
# =============================================================================
if (isTRUE(DO_HUB)) {
  cat("\n[4] Hub-gene tables\n")
  if (is.null(datExpr) || is.null(MEs) || is.null(moduleColors)) {
    message("  Skipped: datExpr/MEs/moduleColors not found in the saved results.")
  } else {
    # Older result bundles saved numeric names despite color-based hub lookup.
    if (all(grepl("^ME[0-9]+$", colnames(MEs)))) {
      colnames(MEs) <- paste0("ME", labels2colors(as.integer(sub("^ME", "", colnames(MEs)))))
    }
    mods <- unique(moduleColors[moduleColors != "grey"])
    if (!is.null(HUB_MODULES)) mods <- intersect(mods, HUB_MODULES)
    if (length(mods) == 0) {
      message("  Skipped: no modules matched HUB_MODULES.")
    } else {
      hdir <- file.path(OUTDIR, "Hub_genes"); dir.create(hdir, showWarnings = FALSE, recursive = TRUE)
      hub_list <- list()
      gene_module <- data.frame(gene = colnames(datExpr), module = moduleColors)
      for (mod in mods) {
        mod_genes <- gene_module$gene[gene_module$module == mod]
        ME <- MEs[, paste0("ME", mod), drop = FALSE]
        kME <- stats::cor(datExpr[, mod_genes, drop = FALSE], ME, use = "pairwise.complete.obs")
        hub <- data.frame(gene = mod_genes, module = mod, kME = as.numeric(kME)) %>%
          dplyr::arrange(dplyr::desc(abs(kME)))
        if (!is.null(HUB_TOP_N)) hub <- utils::head(hub, HUB_TOP_N)
        hub_list[[mod]] <- hub
        write.csv(hub, file.path(hdir, paste0("Hub_genes_", mod, ".csv")), row.names = FALSE)
      }
      write.csv(dplyr::bind_rows(hub_list), file.path(hdir, "Hub_genes_all_modules.csv"), row.names = FALSE)
      cat("  ", length(mods), "module table(s) ->", hdir, "\n")
    }
  }
}

cat("\n========================================\n")
cat("Targeted visualization complete. Outputs under:", OUTDIR, "\n")
cat("========================================\n")
