#!/usr/bin/env Rscript
# Validate the TME deconvolution demo core pipeline.
# Raw counts + gene lengths -> TPM -> ESTIMATE / IOBR deconvolution / ssGSEA.
# Human gene symbols avoid biomaRt, but some IOBR methods may still download
# reference data unless it is already present in the local package cache.
# Run from repository root: Rscript examples/demo_RNAseq_TME/run_demo.R

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_file) == 0 || this_file == "") this_file <- "examples/demo_RNAseq_TME/run_demo.R"
setwd(dirname(this_file))

options(stringsAsFactors = FALSE)

# ---- Parameters ----
RAW_COUNTS_FILE <- "./0-Data/counts.tsv"
RAW_COUNTS_FORMAT <- "tsv"
GENE_COLUMN <- "gene_name"
GENE_LENGTH_COLUMN <- NULL
GENE_LENGTH_UNIT <- "bp"
GENE_START_COL <- "gene_start"
GENE_END_COL <- "gene_end"

META_FILE <- "./0-Data/metadata.csv"
SAMPLE_COLUMN <- "sample"
GROUP_COLUMN <- "condition"
GROUP_LEVELS <- c("Control", "Treatment")

SPECIES <- "human"
GROUP_COLORS <- NULL

RUN_ESTIMATE <- TRUE
RUN_IOBR <- TRUE
IOBR_METHODS <- c("estimate", "cibersort", "epic", "xcell")
IOBR_PERM <- 1000
IOBR_ARRAYS <- FALSE

OUTDIR <- "RNAseq_TME_Deconvolution_Output"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ---- Setup ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(ggpubr)
  library(GSVA)
  library(limma)
})

repo_root <- tryCatch(rprojroot::find_root(rprojroot::is_git_root), error = function(e) getwd())
LIB_DIR <- if (dir.exists("RNAseq_lib")) "RNAseq_lib" else "../../RNAseq_lib"
source(file.path(LIB_DIR, "plot_utils.R"))
source(file.path(LIB_DIR, "io_utils.R"))
source(file.path(LIB_DIR, "tme_utils.R"))
source(file.path(LIB_DIR, "data_utils.R"))
theme_set(theme_publication())

# ---- Load metadata + build TPM expression ----
meta <- read_metadata(
  META_FILE,
  sample_column = SAMPLE_COLUMN,
  required_columns = GROUP_COLUMN,
  group_column = GROUP_COLUMN,
  group_levels = GROUP_LEVELS
)

raw_annot <- read_count_table(RAW_COUNTS_FILE, input_format = RAW_COUNTS_FORMAT)
count_col_names <- detect_count_columns(raw_annot, GENE_COLUMN, annotation_cols = NULL)
counts_mat <- as.matrix(raw_annot[, count_col_names, drop = FALSE])
mode(counts_mat) <- "numeric"
rownames(counts_mat) <- as.character(raw_annot[[GENE_COLUMN]])
validate_count_matrix(counts_mat)

gene_lengths_kb <- extract_gene_lengths(
  raw_annot,
  id_col = GENE_COLUMN,
  length_col = GENE_LENGTH_COLUMN,
  start_col = GENE_START_COL,
  end_col = GENE_END_COL,
  length_unit = GENE_LENGTH_UNIT
)
expr_tpm <- counts_to_tpm(counts_mat, gene_lengths_kb)

common_samples <- intersect(colnames(expr_tpm), meta[[SAMPLE_COLUMN]])
expr_tpm <- expr_tpm[, common_samples, drop = FALSE]
meta <- meta[match(common_samples, meta[[SAMPLE_COLUMN]]), ]

validate_expression_contract(expr_tpm, expected = "tpm")
write.csv(expr_tpm, file.path(OUTDIR, "TPM_matrix.csv"))
cat("Expression:", nrow(expr_tpm), "genes x", ncol(expr_tpm), "samples\n")

# Human symbols: prepare (dedup) the TME matrix without any network conversion.
expr_tme <- prepare_tme_expression(as.data.frame(expr_tpm, check.names = FALSE),
                                   is_log = FALSE, species = "human", verbose = TRUE)
validate_tme_input(expr_tme)
cat("TME expression matrix:", nrow(expr_tme), "genes x", ncol(expr_tme), "samples\n")

if (is.null(GROUP_COLORS) || !all(GROUP_LEVELS %in% names(GROUP_COLORS))) {
  group_colors <- make_group_colors(GROUP_LEVELS)
} else {
  group_colors <- GROUP_COLORS[GROUP_LEVELS]
}
group_df_for_plot <- meta[, c(SAMPLE_COLUMN, GROUP_COLUMN), drop = FALSE]

# ---- Native ESTIMATE ----
# ESTIMATE's filterCommonGenes()/estimateScore() rely on the package's lazy-data
# `common_genes` object. requireNamespace() alone does not resolve that promise,
# so explicitly load it into the environment before calling the ESTIMATE API.
if (RUN_ESTIMATE && requireNamespace("estimate", quietly = TRUE)) {
  utils::data("common_genes", package = "estimate", envir = environment())
  utils::data("SI_geneset", package = "estimate", envir = environment())
  estimate_df <- data.frame(NAME = rownames(expr_tme), Description = NA, expr_tme, check.names = FALSE)
  write.table(estimate_df, file.path(OUTDIR, "estimate_input.gct"), sep = "\t", quote = FALSE, row.names = FALSE)
  estimate::filterCommonGenes(input.f = file.path(OUTDIR, "estimate_input.gct"),
                              output.f = file.path(OUTDIR, "estimate_common_genes.gct"), id = "GeneSymbol")
  estimate::estimateScore(file.path(OUTDIR, "estimate_common_genes.gct"),
                          file.path(OUTDIR, "estimate_scores.gct"), platform = "illumina")
  estimate_scores <- read.table(file.path(OUTDIR, "estimate_scores.gct"), skip = 2, header = TRUE, sep = "\t", check.names = FALSE)
  rownames(estimate_scores) <- estimate_scores$NAME
  estimate_scores <- as.data.frame(t(estimate_scores[, -(1:2)]))
  estimate_scores[[SAMPLE_COLUMN]] <- rownames(estimate_scores)
  write.csv(estimate_scores, file.path(OUTDIR, "ESTIMATE_scores.csv"), row.names = FALSE)
  cat("ESTIMATE scores saved.\n")
}

# ---- IOBR multi-algorithm deconvolution ----
iobr_results <- list()
if (RUN_IOBR && requireNamespace("IOBR", quietly = TRUE)) {
  iobr_results <- run_iobr_deconvolution(
    expr_tme,
    methods = IOBR_METHODS,
    perm = IOBR_PERM,
    arrays = IOBR_ARRAYS,
    id_column = SAMPLE_COLUMN
  )
  for (method in names(iobr_results)) {
    write.csv(iobr_results[[method]], file.path(OUTDIR, paste0("IOBR_", method, ".csv")), row.names = FALSE)
  }
  if (length(iobr_results) >= 2) {
    tme_combined <- combine_tme_results(iobr_results, id_column = SAMPLE_COLUMN)
    write.csv(tme_combined, file.path(OUTDIR, "IOBR_TME_combined.csv"), row.names = FALSE)
  }
  cat("IOBR deconvolution complete:", paste(names(iobr_results), collapse = ", "), "\n")
}

# ---- TME visualization (IOBR) ----
if (RUN_IOBR && "estimate" %in% names(iobr_results)) {
  est_long <- melt_estimate_scores(iobr_results[["estimate"]], id_column = SAMPLE_COLUMN,
                                   group_df = group_df_for_plot,
                                   sample_col = SAMPLE_COLUMN, group_col = GROUP_COLUMN)
  plot_estimate_boxplot_pdf(est_long, group_col = GROUP_COLUMN,
    filename = file.path(OUTDIR, "IOBR_ESTIMATE_scores_boxplot.pdf"),
    title = "ESTIMATE Scores by Group", group_colors = group_colors,
    save_individual = TRUE, individual_prefix = file.path(OUTDIR, "IOBR_ESTIMATE"))
  plot_tme_heatmap_pdf(iobr_results[["estimate"]], meta, group_col = GROUP_COLUMN,
    sample_col = SAMPLE_COLUMN, group_colors = group_colors,
    filename = file.path(OUTDIR, "IOBR_ESTIMATE_heatmap.pdf"), title = "IOBR ESTIMATE Scores")
}

for (method_name in intersect(c("cibersort", "epic"), names(iobr_results))) {
  method_label <- toupper(method_name)
  method_long <- melt_tme_results(iobr_results[[method_name]], id_column = SAMPLE_COLUMN,
    group_df = group_df_for_plot, sample_col = SAMPLE_COLUMN, group_col = GROUP_COLUMN) |>
    dplyr::filter(!grepl("P-value|Correlation|RMSE", .data$cell_type))
  bar_size <- calc_tme_barplot_size(length(unique(method_long[[SAMPLE_COLUMN]])), length(unique(method_long$cell_type)))
  box_size <- calc_tme_boxplot_size(length(unique(method_long$cell_type)))
  prefix <- file.path(OUTDIR, paste0("IOBR_", method_label))
  plot_tme_barplot_pdf(method_long, group_col = GROUP_COLUMN, sample_col = SAMPLE_COLUMN,
    filename = paste0(prefix, "_barplot.pdf"), title = paste(method_label, "Cell Fractions"),
    width = bar_size["width"], height = bar_size["height"])
  plot_tme_boxplot_pdf(method_long, group_col = GROUP_COLUMN, value_col = "fraction",
    filename = paste0(prefix, "_boxplot.pdf"), title = paste(method_label, "Cell Fractions by Group"),
    width = box_size["width"], height = box_size["height"], group_colors = group_colors)
}

if (RUN_IOBR && "xcell" %in% names(iobr_results)) {
  plot_tme_heatmap_pdf(iobr_results[["xcell"]], meta, group_col = GROUP_COLUMN,
    sample_col = SAMPLE_COLUMN, group_colors = group_colors,
    filename = file.path(OUTDIR, "IOBR_xCell_heatmap.pdf"), title = "xCell Scores", width = 10, height = 12)
}

# ---- ssGSEA immune signatures ----
expr_for_ssgsea <- as.data.frame(expr_tme, check.names = FALSE)
row_upper <- toupper(rownames(expr_for_ssgsea))
upper_to_real <- setNames(rownames(expr_for_ssgsea), row_upper)
gs <- lapply(immune_gene_sets, function(x) unique(upper_to_real[intersect(toupper(x), names(upper_to_real))]))
gs <- gs[lengths(gs) >= 2]
stopifnot(length(gs) > 0)

params <- gsvaParam(as.matrix(expr_for_ssgsea), gs, kcdf = "Gaussian", minSize = 2, maxSize = Inf)
ssgsea_scores <- gsva(params, verbose = FALSE)
write.csv(ssgsea_scores, file.path(OUTDIR, "ssGSEA_immune_scores.csv"))

score_df <- as.data.frame(t(ssgsea_scores)) %>% rownames_to_column(SAMPLE_COLUMN) %>% left_join(meta, by = SAMPLE_COLUMN)
score_long <- score_df %>% pivot_longer(cols = names(gs), names_to = "signature", values_to = "score")
score_long[[GROUP_COLUMN]] <- factor(score_long[[GROUP_COLUMN]], levels = GROUP_LEVELS)
score_comparisons <- if (length(GROUP_LEVELS) >= 2) combn(GROUP_LEVELS, 2, simplify = FALSE) else list()
plot_group_boxplot_pdf(
  score_long,
  value_col = "score", group_col = GROUP_COLUMN, facet_col = "signature",
  comparisons = score_comparisons, method = "t.test",
  title = "ssGSEA Immune Signature Scores", ylab = "ssGSEA score",
  group_colors = group_colors,
  filename = file.path(OUTDIR, "ssGSEA_group_boxplot.pdf"), width = 12, height = 8
)

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))

# ---- Assertions ----
expected_files <- c(
  file.path(OUTDIR, "TPM_matrix.csv"),
  file.path(OUTDIR, "ssGSEA_immune_scores.csv"),
  file.path(OUTDIR, "sessionInfo.txt")
)
if (requireNamespace("estimate", quietly = TRUE)) {
  expected_files <- c(expected_files, file.path(OUTDIR, "ESTIMATE_scores.csv"))
}
if (length(iobr_results) > 0) {
  expected_files <- c(expected_files,
    file.path(OUTDIR, "IOBR_TME_combined.csv"),
    file.path(OUTDIR, "IOBR_ESTIMATE_heatmap.pdf"))
}
missing_files <- expected_files[!file.exists(expected_files)]
if (length(missing_files) > 0) {
  stop("TME demo FAILED. Missing expected outputs:\n", paste(missing_files, collapse = "\n"))
}
stopifnot(nrow(ssgsea_scores) > 0, ncol(ssgsea_scores) == length(common_samples))

cat("\n========================================\n")
cat("TME deconvolution demo PASSED.\n")
cat("ssGSEA signatures scored:", nrow(ssgsea_scores), "\n")
cat("Outputs saved to", OUTDIR, "\n")
cat("========================================\n")
