#!/usr/bin/env Rscript
# Validate the TimeCourse demo core pipeline.
# VST expression across time points -> Mfuzz soft clustering -> per-cluster ORA,
# plus raw-count time-point-vs-baseline DESeq2.
# Run from repository root: Rscript examples/demo_RNAseq_TimeCourse/run_demo.R

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_file) == 0 || this_file == "") this_file <- "examples/demo_RNAseq_TimeCourse/run_demo.R"
setwd(dirname(this_file))

options(stringsAsFactors = FALSE)

# ---- Parameters ----
EXPR_FILE <- "./0-Data/vsd_matrix.csv"
META_FILE <- "./0-Data/colData.csv"
GENE_COLUMN <- "gene_name"
SAMPLE_COLUMN <- "sample"
TIME_COLUMN <- "time"
GROUP_COLUMN <- "condition"
TIME_LEVELS <- c("Day0", "Day7", "Day14", "Day21")

RUN_MFUZZ <- TRUE
MFUZZ_N_CLUSTERS <- 4
MFUZZ_MIN_ACORE <- 0.7
MFUZZ_SEED <- 2025

RAW_COUNTS_FILE <- "./0-Data/raw_counts.tsv"
COUNT_META_FILE <- "./0-Data/metadata.csv"
COUNT_GENE_COL <- "gene_name"
COUNT_SAMPLE_COL <- "sample"
COUNT_BIOTYPE_COL <- NULL
COUNT_BIOTYPE_FILTER <- "protein_coding"

SUBJECT_COL <- NULL
RUN_TIMEPOINT_DEG <- TRUE
BASELINE_TIME <- "Day0"
DEG_PADJ_CUTOFF <- 0.05
DEG_LFC_CUTOFF <- 0.5
MIN_COUNT <- 10

OUTDIR <- "RNAseq_TimeCourse_Output"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ---- Setup ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(clusterProfiler)
  library(DESeq2)
  library(org.Hs.eg.db)
})
# Mfuzz is optional for this demo; load it only when available so the rest of
# the pipeline (aggregation, time-point DEG) still runs on minimal systems.
if (requireNamespace("Mfuzz", quietly = TRUE)) {
  suppressPackageStartupMessages(library(Mfuzz))
}

LIB_DIR <- if (dir.exists("RNAseq_lib")) "RNAseq_lib" else "../../RNAseq_lib"
source(file.path(LIB_DIR, "plot_utils.R"))
source(file.path(LIB_DIR, "io_utils.R"))
source(file.path(LIB_DIR, "data_utils.R"))
source(file.path(LIB_DIR, "deg_utils.R"))
source(file.path(LIB_DIR, "enrichment_utils.R"))
source(file.path(LIB_DIR, "timecourse_utils.R"))
theme_set(theme_publication())

# ---- Load expression + metadata ----
expr <- read_expression_matrix(EXPR_FILE, gene_column = GENE_COLUMN)
meta <- read_metadata(
  META_FILE,
  sample_column = SAMPLE_COLUMN,
  required_columns = c(SAMPLE_COLUMN, TIME_COLUMN, GROUP_COLUMN),
  time_column = TIME_COLUMN,
  time_levels = TIME_LEVELS,
  group_column = GROUP_COLUMN
)
validate_samples_match(colnames(expr), meta[[SAMPLE_COLUMN]], strict_order = TRUE)
expr <- expr[, meta[[SAMPLE_COLUMN]], drop = FALSE]
validate_expression_contract(expr, expected = "vst")
cat("Expression:", nrow(expr), "genes x", ncol(expr), "samples\n")

# ---- Aggregate to per-time-point means ----
expr_mean <- aggregate_expr_by_group(expr, meta[[TIME_COLUMN]])
expr_mean <- expr_mean[, TIME_LEVELS, drop = FALSE]
write.csv(expr_mean, file.path(OUTDIR, "mean_expression_by_time.csv"))

# ---- Mfuzz clustering ----
mfuzz_available <- RUN_MFUZZ && requireNamespace("Mfuzz", quietly = TRUE)
if (RUN_MFUZZ && !mfuzz_available) {
  message("Package 'Mfuzz' is not available; skipping Mfuzz clustering steps.")
}
cluster_df <- NULL
if (mfuzz_available) {
  eset <- prepare_mfuzz_eset(expr_mean)
  mfuzz_result <- run_mfuzz(eset, n_clusters = MFUZZ_N_CLUSTERS, seed = MFUZZ_SEED)
  cluster_df <- extract_mfuzz_clusters(mfuzz_result, eset = eset, min_acore = MFUZZ_MIN_ACORE)
  write_mfuzz_cluster_table(cluster_df, file.path(OUTDIR, "mfuzz_clusters.csv"))
  cat("Mfuzz cluster sizes:\n")
  print(summarize_mfuzz_clusters(cluster_df))

  plot_mfuzz_trends_pdf(
    eset, mfuzz_result,
    filename = file.path(OUTDIR, "mfuzz_trends.pdf"),
    time_labels = TIME_LEVELS, width = 12, height = 9
  )

  core_df <- cluster_df[cluster_df$core_gene, ]
  if (nrow(core_df) >= 2) {
    group_colors <- make_group_colors(TIME_LEVELS)
    plot_timecourse_heatmap_pdf(
      expr, core_df,
      group_vec = meta[[TIME_COLUMN]],
      group_levels = TIME_LEVELS,
      group_colors = group_colors,
      filename = file.path(OUTDIR, "mfuzz_core_heatmap.pdf"),
      width = 9, height = 12
    )
  }

  # ---- ORA per cluster ----
  universe <- map_symbols_to_entrez(rownames(expr), org.Hs.eg.db)$ENTREZID
  ora_results <- run_mfuzz_cluster_ora(cluster_df, org_db = org.Hs.eg.db, universe = universe)
  for (cl_name in names(ora_results)) {
    prefix <- file.path(OUTDIR, paste0("GO_ORA_", cl_name))
    write.csv(as.data.frame(ora_results[[cl_name]]), paste0(prefix, ".csv"), row.names = FALSE)
    plot_enrich_suite_pdf(ora_results[[cl_name]], prefix, cl_name)
  }
}

# ---- Time-point vs baseline DESeq2 ----
if (RUN_TIMEPOINT_DEG) {
  rawcount <- read_count_table(RAW_COUNTS_FILE, "tsv")
  count_col_names <- detect_count_columns(rawcount, COUNT_GENE_COL, NULL)
  count_meta <- read.csv(COUNT_META_FILE, check.names = FALSE)
  count_samples <- intersect(count_col_names, count_meta[[COUNT_SAMPLE_COL]])
  rawcount <- rawcount[, c(COUNT_GENE_COL, count_samples), drop = FALSE]
  count_meta <- count_meta[match(count_samples, count_meta[[COUNT_SAMPLE_COL]]), ]

  countData <- build_count_matrix(
    rawcount, COUNT_GENE_COL, count_samples, count_meta[[COUNT_SAMPLE_COL]],
    duplicate_report_file = file.path(OUTDIR, "1-DEG_Timepoint", "Duplicated_gene_symbols.csv")
  )
  countData <- filter_low_count_genes(countData, count_meta[[TIME_COLUMN]], MIN_COUNT)$count_data

  col_data <- data.frame(sample = count_meta[[COUNT_SAMPLE_COL]], stringsAsFactors = FALSE)
  col_data[[TIME_COLUMN]] <- count_meta[[TIME_COLUMN]]
  # Only include the group/condition column in the design when it actually
  # varies across samples; a single-level condition breaks the DESeq2 design.
  condition_col <- NULL
  if (GROUP_COLUMN %in% colnames(count_meta) &&
      length(unique(count_meta[[GROUP_COLUMN]])) > 1) {
    col_data[[GROUP_COLUMN]] <- count_meta[[GROUP_COLUMN]]
    condition_col <- GROUP_COLUMN
  }
  rownames(col_data) <- col_data$sample

  tp_res_list <- run_timepoint_vs_baseline_deseq2(
    count_data = countData,
    col_data = col_data,
    time_col = TIME_COLUMN,
    baseline_time = BASELINE_TIME,
    condition_col = condition_col,
    subject_col = SUBJECT_COL,
    alpha = max(DEG_PADJ_CUTOFF, 0.05)
  )

  tp_deg_dir <- file.path(OUTDIR, "1-DEG_Timepoint")
  tp_summary <- write_timepoint_deg_results(
    tp_res_list, outdir = tp_deg_dir,
    pvalue_column = "padj", lfc_column = "log2FoldChange_shrunken"
  )

  plot_timepoint_deg_summary_pdf(
    tp_summary, filename = file.path(OUTDIR, "3-Visualization", "Timepoint_DEG_summary.pdf")
  )

  for (comp_name in names(tp_res_list)) {
    plot_volcano_pdf(
      tp_res_list[[comp_name]], comp_name = comp_name,
      pvalue_thresh = DEG_PADJ_CUTOFF, log2fc_thresh = DEG_LFC_CUTOFF,
      pvalue_column = "padj", lfc_column = "log2FoldChange_shrunken",
      filename = file.path(OUTDIR, "3-Visualization", paste0("Volcano_", comp_name, ".pdf"))
    )
  }
  cat("Time-point vs baseline DEG complete.\n")
}

save.image(file = file.path(OUTDIR, "timecourse_workspace.Rdata"))
writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))

# ---- Assertions ----
# Core outputs that must always be produced regardless of optional packages.
expected_files <- c(
  file.path(OUTDIR, "mean_expression_by_time.csv"),
  file.path(OUTDIR, "sessionInfo.txt")
)
if (mfuzz_available) {
  expected_files <- c(expected_files, file.path(OUTDIR, "mfuzz_clusters.csv"))
}
if (RUN_TIMEPOINT_DEG) {
  expected_files <- c(expected_files,
    file.path(OUTDIR, "1-DEG_Timepoint", "Timepoint_DEG_summary.csv"))
}
missing_files <- expected_files[!file.exists(expected_files)]
if (length(missing_files) > 0) {
  stop("TimeCourse demo FAILED. Missing expected outputs:\n", paste(missing_files, collapse = "\n"))
}
if (mfuzz_available) {
  stopifnot(nrow(cluster_df) > 0, length(unique(cluster_df$cluster)) >= 2)
}
if (RUN_TIMEPOINT_DEG && isTRUE(capabilities("cairo"))) {
  deg_pdf <- file.path(OUTDIR, "3-Visualization", "Timepoint_DEG_summary.pdf")
  if (!file.exists(deg_pdf)) warning("TimeCourse demo: expected PDF not generated: ", deg_pdf)
}

cat("\n========================================\n")
cat("TimeCourse demo PASSED.\n")
if (mfuzz_available) {
  cat("Mfuzz clusters:", length(unique(cluster_df$cluster)), "| core genes:", sum(cluster_df$core_gene), "\n")
} else {
  cat("Mfuzz clustering skipped (package unavailable); time-point DEG ran.\n")
}
cat("Outputs saved to", OUTDIR, "\n")
cat("========================================\n")
