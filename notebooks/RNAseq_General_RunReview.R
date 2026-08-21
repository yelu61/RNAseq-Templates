# ---
# jupyter:
#   jupytext:
#     formats: ipynb,R:percent
#     text_representation:
#       extension: .R
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---

# %% [markdown]
# # RNAseq General — existing-run review
#
# This notebook is the read-only companion to the complete analysis notebook.
# It reviews one completed CLI run, replots selected genes, and optionally
# recalculates user-supplied gene-set scores. It never fits the core statistical
# model and never writes into the source run directory.

# %% [markdown]
# ## 1. Review parameters — edit here

# %%
PROJECT_ROOT <- Sys.getenv("RNASEQ_PROJECT_ROOT", unset = normalizePath(".", mustWork = TRUE))
RUN_DIR <- Sys.getenv("RNASEQ_RUN_DIR",
                      unset = file.path(PROJECT_ROOT, "analysis", "runs", "<run_id>"))
REVIEW_OUTDIR <- Sys.getenv("RNASEQ_REVIEW_OUTDIR",
                            unset = file.path(PROJECT_ROOT, "analysis", "notebook_output", "run_review", "<run_id>"))

PRIMARY_THRESHOLD <- "standard"
COMPARISON <- NULL                 # NULL = first comparison found in the run
KEY_GENES <- c("Tnf", "Il1b", "Cxcl10", "Mrc1", "Arg1")

RUN_CUSTOM_GSVA <- FALSE
CUSTOM_GENE_SETS <- list(
  Inflammatory_response = c("Il1b", "Il6", "Tnf", "Cxcl1", "Cxcl2", "Ccl2"),
  Anti_inflammatory = c("Il10", "Tgfb1", "Mrc1", "Arg1", "Cd163")
)

# %% [markdown]
# ## 2. Read-only guard and shared plotting library

# %%
run_dir <- normalizePath(path.expand(RUN_DIR), mustWork = TRUE)
review_outdir <- normalizePath(path.expand(REVIEW_OUTDIR), mustWork = FALSE)
if (identical(run_dir, review_outdir) || startsWith(review_outdir, paste0(run_dir, .Platform$file.sep))) {
  stop("REVIEW_OUTDIR must be outside RUN_DIR; the completed run is read-only.")
}
required <- c("1-DEG/vsd_matrix.csv", "1-DEG/colData.csv", "Analysis_summary.txt")
missing <- required[!file.exists(file.path(run_dir, required))]
if (length(missing)) stop("RUN_DIR is incomplete; missing: ", paste(missing, collapse = ", "))
dir.create(review_outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(review_outdir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(review_outdir, "tables"), recursive = TRUE, showWarnings = FALSE)

lib_dir <- Sys.getenv("RNASEQ_LIB_DIR", unset = file.path(PROJECT_ROOT, "RNAseq_lib"))
if (!dir.exists(lib_dir)) stop("Set RNASEQ_LIB_DIR to the RNAseq-Templates/RNAseq_lib directory.")
source(file.path(lib_dir, "plot_utils.R"))
source(file.path(lib_dir, "pathway_utils.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
theme_set(theme_publication())

# %% [markdown]
# ## 3. Load the saved expression matrix and design

# %%
vsd_df <- read.csv(file.path(run_dir, "1-DEG", "vsd_matrix.csv"),
                   check.names = FALSE, stringsAsFactors = FALSE)
gene_col <- colnames(vsd_df)[[1]]
vsd_mat <- as.matrix(vsd_df[, -1, drop = FALSE])
storage.mode(vsd_mat) <- "double"
rownames(vsd_mat) <- vsd_df[[gene_col]]

coldata <- read.csv(file.path(run_dir, "1-DEG", "colData.csv"),
                    check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(coldata) != ncol(vsd_mat)) stop("colData rows do not match VST samples.")
if (!"sample" %in% colnames(coldata)) coldata$sample <- colnames(vsd_mat)
if (!"condition" %in% colnames(coldata)) stop("colData.csv has no condition column.")
coldata <- coldata[match(colnames(vsd_mat), coldata$sample), , drop = FALSE]
group_levels <- unique(as.character(coldata$condition))
group_colors <- make_group_colors(group_levels)

manifest_file <- file.path(run_dir, "run_manifest.csv")
run_manifest <- if (file.exists(manifest_file)) read.csv(manifest_file, stringsAsFactors = FALSE) else NULL
list(run = basename(run_dir), samples = ncol(vsd_mat), genes = nrow(vsd_mat),
     groups = table(coldata$condition), manifest = run_manifest)

# %% [markdown]
# ## 4. Review saved differential-expression results

# %%
deg_files <- list.files(file.path(run_dir, "1-DEG", PRIMARY_THRESHOLD),
                        pattern = "^DEG_results_.*[.]csv$", full.names = TRUE)
if (!length(deg_files)) stop("No DEG result files found for threshold: ", PRIMARY_THRESHOLD)
comparison_names <- sub("[.]csv$", "", sub("^DEG_results_", "", basename(deg_files)))
if (is.null(COMPARISON)) COMPARISON <- comparison_names[[1]]
deg_file <- deg_files[match(COMPARISON, comparison_names)]
if (is.na(deg_file)) stop("Unknown COMPARISON. Available: ", paste(comparison_names, collapse = ", "))
deg <- read.csv(deg_file, stringsAsFactors = FALSE, check.names = FALSE)
deg_counts <- table(factor(deg$significance, levels = c("Up", "Down", "Not_Sig")))
top_deg <- deg[order(deg$padj, -abs(deg$log2FoldChange_raw)), , drop = FALSE]
head(top_deg, 20)
deg_counts

# %% [markdown]
# ## 5. Selected-gene visualization from the saved VST matrix

# %%
genes_present <- intersect(KEY_GENES, rownames(vsd_mat))
if (length(genes_present)) {
  plot_gene_expression_pdf(
    vsd_mat, genes = genes_present,
    group = setNames(as.character(coldata$condition), coldata$sample),
    group_levels = group_levels, group_colors = group_colors,
    filename = file.path(review_outdir, "figures", "Selected_genes.pdf"),
    plot = "violin", facet = TRUE,
    comparisons = if (length(group_levels) > 1) combn(group_levels, 2, simplify = FALSE) else list(),
    method = "t.test", p_adjust_method = "BH",
    title = paste("Selected genes —", COMPARISON), ylab = "VST expression"
  )
}
data.frame(requested = KEY_GENES, present = KEY_GENES %in% rownames(vsd_mat))

# %% [markdown]
# ## 6. Inspect existing pathway and TME outputs

# %%
read_if_present <- function(...) {
  path <- file.path(run_dir, ...)
  if (file.exists(path)) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE) else NULL
}
gsva_saved <- read_if_present("2-GSEA", "GSVA_group_comparison.csv")
gsea_go <- read_if_present("2-GSEA", paste0("GSEA_GO_", COMPARISON, ".csv"))
tme_files <- if (dir.exists(file.path(run_dir, "4-TME"))) {
  list.files(file.path(run_dir, "4-TME"), recursive = TRUE)
} else character(0)

if (!is.null(gsva_saved)) head(gsva_saved[order(gsva_saved$p.adj), ], 20)
if (!is.null(gsea_go)) head(gsea_go[order(gsea_go$p.adjust), ], 20)
tme_files

# %% [markdown]
# ## 7. Optional custom GSVA derived from the saved VST matrix
#
# This is a new exploratory derivative. It is written only to REVIEW_OUTDIR and
# must not be confused with the custom gene sets frozen in the source run.

# %%
if (isTRUE(RUN_CUSTOM_GSVA)) {
  custom_scores <- score_gene_sets(vsd_mat, CUSTOM_GENE_SETS, method = "gsva", min_size = 3)
  custom_stats <- pathway_group_comparison(
    custom_scores,
    group = setNames(as.character(coldata$condition), coldata$sample),
    group_levels = group_levels,
    comparisons = if (length(group_levels) > 1) combn(group_levels, 2, simplify = FALSE) else list(),
    method = "t.test", p_adjust_method = "BH"
  )
  write.csv(data.frame(gene_set = rownames(custom_scores), custom_scores, check.names = FALSE),
            file.path(review_outdir, "tables", "custom_GSVA_scores.csv"), row.names = FALSE)
  write.csv(custom_stats, file.path(review_outdir, "tables", "custom_GSVA_group_comparison.csv"),
            row.names = FALSE)
  plot_expression_heatmap_pdf(
    custom_scores, file.path(review_outdir, "figures", "custom_GSVA_heatmap.pdf"),
    title = "Review-only custom GSVA", group = as.character(coldata$condition),
    group_levels = group_levels, group_colors = group_colors,
    scale_rows = FALSE, show_row_names = TRUE, show_column_names = TRUE,
    width = 8, height = max(4, 0.6 * nrow(custom_scores) + 2)
  )
  custom_stats
}

# %% [markdown]
# ## 8. Review provenance

# %%
provenance <- data.frame(
  source_run_id = basename(run_dir),
  source_manifest_md5 = if (file.exists(manifest_file)) unname(tools::md5sum(manifest_file)) else NA_character_,
  comparison = COMPARISON,
  threshold = PRIMARY_THRESHOLD,
  custom_gsva_ran = isTRUE(RUN_CUSTOM_GSVA),
  reviewed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  stringsAsFactors = FALSE
)
write.csv(provenance, file.path(review_outdir, "review_provenance.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(review_outdir, "sessionInfo.txt"))
provenance
