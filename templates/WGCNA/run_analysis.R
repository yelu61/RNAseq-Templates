#!/usr/bin/env Rscript
# =============================================================================
# RNAseq_WGCNA — non-interactive analysis pipeline
# =============================================================================
# Runs the SAME workflow as notebooks/RNAseq_WGCNA_Template.ipynb, driven by a
# config file instead of a parameter cell. Intended use:
#
#   1. Copy templates/WGCNA/ (config.R + run_analysis.R) into your project.
#   2. Edit config.R (expression + trait paths, filtering, network params).
#   3. Run:  Rscript run_analysis.R
#      or:   Rscript run_analysis.R path/to/other_config.R
#
# LIB_DIR resolution: if ./RNAseq_lib exists next to this script it is used;
# otherwise the repository root is located via rprojroot and its RNAseq_lib is
# used. Set the environment variable RNASEQ_LIB_DIR to override explicitly.
#
# There is no wgcna_utils module — this runner calls the WGCNA package directly
# and sources only plot_utils.R, data_utils.R and report_utils.R from RNAseq_lib.
# =============================================================================

options(stringsAsFactors = FALSE)

# ---- Locate this script and the config file ----------------------------------
invocation_dir <- normalizePath(getwd(), mustWork = TRUE)
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", cmd_args, value = TRUE))
file_arg <- gsub("~\\+~", " ", file_arg)
script_dir <- if (length(file_arg) > 0 && nzchar(file_arg)) dirname(normalizePath(file_arg[1])) else invocation_dir

user_args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(user_args) >= 1) normalizePath(path.expand(user_args[1]), mustWork = FALSE) else file.path(script_dir, "config.R")
if (!file.exists(config_path)) {
  stop("Config file not found: ", config_path,
       "\nProvide one as: Rscript run_analysis.R path/to/config.R")
}
config_path <- normalizePath(config_path, mustWork = TRUE)

cat("========================================\n")
cat("RNAseq_WGCNA — run_analysis.R\n")
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
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  stop("Package 'WGCNA' is required. Install with BiocManager::install('WGCNA') ",
       "or install.packages('WGCNA').")
}
suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
})
WGCNA::allowWGCNAThreads()

for (f in c("plot_utils.R", "data_utils.R", "report_utils.R")) {
  source(file.path(lib_dir, f))
}
theme_set(theme_publication())

# ---- Output directories + config snapshot ------------------------------------
# WGCNA tables/figures go under OUTDIR (default "5-WGCNA"). The config snapshot,
# run summary and session info are run-root artifacts, written next to OUTDIR
# (the parent of OUTDIR) so the numbered layout stays grouped.
RUN_ROOT   <- dirname(OUTDIR)
CONFIG_DIR <- file.path(RUN_ROOT, "0-Config")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CONFIG_DIR, showWarnings = FALSE, recursive = TRUE)

config_objects <- c(
  "EXPR_FILE", "TRAIT_FILE", "GENE_COLUMN", "SAMPLE_COLUMN", "GROUP_COLUMN",
  "MIN_MAD_QUANTILE", "NETWORK_TYPE", "POWER_VECTOR", "SOFT_POWER",
  "MIN_MODULE_SIZE", "MERGE_CUT_HEIGHT", "TARGET_MODULES", "OUTDIR",
  "GENERATE_HTML_REPORT", "REPORT_TITLE"
)
config_lines <- c(
  "# RNAseq_WGCNA analysis configuration snapshot",
  paste0("# Saved: ", Sys.time()),
  unlist(lapply(config_objects, function(obj) {
    if (exists(obj, inherits = TRUE)) c(paste0("\n", obj, " <- "), capture.output(dput(get(obj, inherits = TRUE)))) else NULL
  }))
)
writeLines(config_lines, file.path(CONFIG_DIR, "analysis_config_used.R"))

# =============================================================================
# 3. Load expression and traits
# =============================================================================
expr <- read_expression_matrix(EXPR_FILE, gene_column = GENE_COLUMN)

traits <- read_metadata(
  TRAIT_FILE,
  sample_column = SAMPLE_COLUMN,
  required_columns = NULL,
  group_column = GROUP_COLUMN
)

validate_samples_match(colnames(expr), traits[[SAMPLE_COLUMN]], strict_order = TRUE)
expr <- expr[, traits[[SAMPLE_COLUMN]], drop = FALSE]
rownames(traits) <- traits[[SAMPLE_COLUMN]]

# WGCNA requires variance-stabilized/normalized expression derived from raw counts.
scale_info <- validate_expression_contract(expr, expected = "vst")

cat("Expression:", nrow(expr), "genes x", ncol(expr), "samples\n")

# =============================================================================
# 4. Gene filtering and sample QC
# =============================================================================
gene_mad <- apply(expr, 1, mad, na.rm = TRUE)
expr <- expr[gene_mad >= quantile(gene_mad, MIN_MAD_QUANTILE, na.rm = TRUE), , drop = FALSE]
datExpr <- t(as.matrix(expr))

gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  traits <- traits[rownames(datExpr), , drop = FALSE]
}

sampleTree <- hclust(dist(datExpr), method = "average")
plot_wgcna_sample_tree_pdf(sampleTree, file.path(OUTDIR, "Sample_clustering.pdf"))
cat("WGCNA input:", nrow(datExpr), "samples x", ncol(datExpr), "genes\n")

# =============================================================================
# 5. Soft-threshold selection
# =============================================================================
sft <- pickSoftThreshold(datExpr, powerVector = POWER_VECTOR, networkType = NETWORK_TYPE, verbose = 0)

if (!is.null(SOFT_POWER)) {
  soft_power <- SOFT_POWER
  cat("Using SOFT_POWER override from config:", soft_power, "\n")
} else {
  soft_power <- sft$powerEstimate
  if (is.na(soft_power)) {
    fit_df <- sft$fitIndices
    soft_power <- fit_df$Power[which.max(fit_df$SFT.R.sq)]
    message("No automatic power estimate; using max scale-free fit power: ", soft_power)
  }
}

plot_wgcna_soft_threshold_pdf(
  sft, selected_power = soft_power,
  filename = file.path(OUTDIR, "Soft_threshold_selection.pdf")
)
cat("Selected soft power:", soft_power, "\n")

# =============================================================================
# 6. Network construction and module detection
# =============================================================================
net <- blockwiseModules(
  datExpr,
  power = soft_power,
  networkType = NETWORK_TYPE,
  TOMType = NETWORK_TYPE,
  minModuleSize = MIN_MODULE_SIZE,
  reassignThreshold = 0,
  mergeCutHeight = MERGE_CUT_HEIGHT,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 0
)
moduleColors <- labels2colors(net$colors)
MEs <- orderMEs(net$MEs)

write.csv(data.frame(gene = colnames(datExpr), module = moduleColors),
          file.path(OUTDIR, "WGCNA_gene_modules.csv"), row.names = FALSE)
saveRDS(list(net = net, moduleColors = moduleColors, MEs = MEs, datExpr = datExpr, traits = traits),
        file.path(OUTDIR, "WGCNA_network.rds"))

plot_wgcna_module_dendrogram_pdf(
  net, moduleColors, file.path(OUTDIR, "Module_dendrogram.pdf")
)
print(table(moduleColors))

# =============================================================================
# 7. Module-trait correlation
# =============================================================================
trait_numeric <- encode_wgcna_traits(traits, sample_column = SAMPLE_COLUMN)
if (ncol(trait_numeric) == 0) {
  stop("No usable trait columns found in TRAIT_FILE after dropping the sample ",
       "identifier mirrors, all-NA and invariant columns. Provide at least one ",
       "varying numeric or categorical ",
       "trait (besides '", SAMPLE_COLUMN, "') to correlate modules against.")
}
moduleTraitCor <- cor(MEs, trait_numeric, use = "p")
moduleTraitP <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

write.csv(moduleTraitCor, file.path(OUTDIR, "Module_trait_correlation.csv"))
write.csv(moduleTraitP, file.path(OUTDIR, "Module_trait_pvalue.csv"))

plot_wgcna_module_trait_heatmap_pdf(
  moduleTraitCor, moduleTraitP, file.path(OUTDIR, "Module_trait_heatmap.pdf")
)

# =============================================================================
# 8. Hub-gene export
# =============================================================================
gene_module <- data.frame(gene = colnames(datExpr), module = moduleColors)
all_modules <- unique(moduleColors[moduleColors != "grey"])
if (!is.null(TARGET_MODULES)) all_modules <- intersect(all_modules, TARGET_MODULES)

hub_list <- list()
for (mod in all_modules) {
  mod_genes <- gene_module$gene[gene_module$module == mod]
  ME <- MEs[[paste0("ME", mod)]]
  kME <- cor(datExpr[, mod_genes, drop = FALSE], ME, use = "p")
  hub <- data.frame(gene = mod_genes, module = mod, kME = as.numeric(kME)) %>%
    dplyr::arrange(dplyr::desc(abs(kME)))
  hub_list[[mod]] <- hub
  write.csv(hub, file.path(OUTDIR, paste0("Hub_genes_", mod, ".csv")), row.names = FALSE)
}
write.csv(dplyr::bind_rows(hub_list), file.path(OUTDIR, "Hub_genes_all_modules.csv"), row.names = FALSE)

# =============================================================================
# 9. Save results, summary, session info
# =============================================================================
# Targeted save so visualize_results.R can re-plot the dendrogram / module-trait
# heatmap / hub tables without re-running blockwiseModules.
save(net, MEs, moduleColors, datExpr, traits, trait_numeric, moduleTraitCor,
     moduleTraitP, sft, soft_power, all_modules, hub_list,
     NETWORK_TYPE, POWER_VECTOR, MIN_MODULE_SIZE, MERGE_CUT_HEIGHT,
     SAMPLE_COLUMN, GROUP_COLUMN,
     file = file.path(OUTDIR, "WGCNA_results.Rdata"))

module_counts <- table(moduleColors)
summary_report <- paste0(
  "========================================\n",
  "WGCNA Analysis Summary\n",
  "========================================\n\n",
  "1. Data Overview\n",
  "   - Expression file: ", EXPR_FILE, "\n",
  "   - Trait file: ", TRAIT_FILE, "\n",
  "   - Samples used: ", nrow(datExpr), "\n",
  "   - Genes after MAD filter (quantile ", MIN_MAD_QUANTILE, "): ", ncol(datExpr), "\n",
  "   - Trait columns modelled: ", ncol(trait_numeric), "\n\n",
  "2. Network\n",
  "   - Network type: ", NETWORK_TYPE, "\n",
  "   - Soft-threshold power: ", soft_power, "\n",
  "   - Min module size: ", MIN_MODULE_SIZE, " | merge cut height: ", MERGE_CUT_HEIGHT, "\n",
  "   - Modules detected (incl. grey): ", length(module_counts), "\n"
)
for (mod in names(module_counts)) {
  summary_report <- paste0(summary_report, "      - ", mod, ": ", module_counts[[mod]], " genes\n")
}
summary_report <- paste0(summary_report, "\n3. Output Files\n",
  "   - ", file.path(CONFIG_DIR, "analysis_config_used.R"), "\n",
  "   - ", OUTDIR, " (WGCNA_network.rds, gene-module + module-trait tables, hub-gene CSVs, PDFs, results Rdata)\n",
  "\n========================================\n")
cat(summary_report)
writeLines(summary_report, file.path(RUN_ROOT, "Analysis_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(RUN_ROOT, "sessionInfo.txt"))

# =============================================================================
# 10. HTML report (optional)
# =============================================================================
if (isTRUE(GENERATE_HTML_REPORT)) {
  has_quarto    <- requireNamespace("quarto", quietly = TRUE) && nzchar(Sys.which("quarto"))
  has_rmarkdown <- requireNamespace("rmarkdown", quietly = TRUE)
  if (!has_quarto && !has_rmarkdown) {
    message("Skipping HTML report: install the quarto CLI (+ R 'quarto' package) or the 'rmarkdown' package.")
  } else {
    report_path <- tryCatch(
      render_analysis_report(outdir = RUN_ROOT, report_file = file.path(RUN_ROOT, "RNAseq_report.html"),
                             params = list(title = REPORT_TITLE, author = Sys.info()[["user"]])),
      error = function(e) { message("HTML report failed: ", conditionMessage(e)); NULL }
    )
    if (!is.null(report_path)) cat("HTML report written to:", report_path, "\n")
  }
}

cat("\n========================================\n")
cat("RNAseq_WGCNA analysis COMPLETE.\n")
cat("Outputs under:", normalizePath(OUTDIR), "\n")
cat("========================================\n")
