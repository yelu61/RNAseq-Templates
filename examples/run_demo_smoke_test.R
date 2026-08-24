#!/usr/bin/env Rscript
# Smoke test for RNAseq-Templates.
# Runs the General helper-level pipeline plus every production runner on bundled
# demo data and checks that all expected output files are produced.
#
# Usage from repository root:
#   Rscript examples/run_demo_smoke_test.R

options(stringsAsFactors = FALSE)

# Locate repository root and helper files
# When run with Rscript, script_path is the first command-line argument.
script_path <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_path[grep("^--file=", script_path)])
if (length(script_path) == 0 || script_path == "") script_path <- "examples/run_demo_smoke_test.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."))
setwd(repo_root)
lib_dir <- file.path(repo_root, "RNAseq_lib")
if (!dir.exists(lib_dir)) {
  stop("RNAseq_lib directory not found at: ", lib_dir)
}

cat("========================================\n")
cat("RNAseq-Templates Demo Smoke Test\n")
cat("Working directory:", getwd(), "\n")
cat("RNAseq_lib:", lib_dir, "\n")
cat("========================================\n\n")

# ---- Check required packages ----
required_pkgs <- c(
  "DESeq2", "clusterProfiler", "org.Mm.eg.db", "ComplexHeatmap",
  "circlize", "matrixStats", "EnhancedVolcano", "GSVA", "msigdbr",
  "ggplot2", "dplyr", "tidyr", "RColorBrewer", "ggrepel", "DOSE",
  "enrichplot", "ggpubr", "tidyverse", "data.table", "pheatmap", "ashr", "BiocParallel"
)
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing required packages. Run install_dependencies.R first:\n",
       paste(missing, collapse = ", "))
}

# ---- Load helpers ----
source(file.path(lib_dir, "plot_utils.R"))
source(file.path(lib_dir, "io_utils.R"))
source(file.path(lib_dir, "deg_utils.R"))
source(file.path(lib_dir, "enrichment_utils.R"))
source(file.path(lib_dir, "design_utils.R"))
source(file.path(lib_dir, "data_utils.R"))

# Load organism annotation package
suppressPackageStartupMessages(library(org.Mm.eg.db))
org_db <- org.Mm.eg.db

# ---- Demo parameters ----
SPECIES <- "mouse"
INPUT_FILE <- file.path(repo_root, "examples/demo_data/demo_counts.tsv")
INPUT_FORMAT <- "tsv"
GENE_NAME_COL <- "gene_name"
BIOTYPE_COL <- "gene_biotype"
BIOTYPE_FILTER <- "protein_coding"
COUNT_COLS <- NULL

SAMPLE_NAMES <- c("Control_1", "Control_2", "Control_3",
                  "Treatment_1", "Treatment_2", "Treatment_3")
GROUPS <- c(rep("Control", 3), rep("Treatment", 3))
GROUP_LEVELS <- c("Control", "Treatment")

COMPARISONS <- list(c("Treatment_vs_Control", "Treatment", "Control"))

THRESHOLD_GRID <- data.frame(
  name = c("strict", "standard", "loose"),
  p_cutoff = c(0.01, 0.05, 0.10),
  log2fc = c(1.5, 1.0, 0.5),
  stringsAsFactors = FALSE
)
DEFAULT_THRESHOLD <- "standard"

MIN_COUNT <- 10
DESIGN_FORMULA <- ~ condition
DEG_PVALUE_COLUMN <- "padj"
DEG_LFC_COLUMN <- "log2FoldChange_raw"
GSEA_RANK_COLUMN <- "stat"
PAIRWISE_TEST_METHOD <- "t.test"
PAIRWISE_P_ADJUST_METHOD <- "BH"

MIN_LIBRARY_SIZE <- NULL
MIN_DETECTED_GENES <- NULL
MAX_ZERO_FRACTION <- NULL
MIN_MEDIAN_CORRELATION_Z <- -3
SAMPLE_EXCLUDE <- character(0)

OUTDIR <- file.path(repo_root, "examples/demo_RNAseq_General")
setwd(OUTDIR)

# Ensure output directories exist
dir.create("0-Config", showWarnings = FALSE, recursive = TRUE)
dir.create("1-DEG", showWarnings = FALSE, recursive = TRUE)
dir.create("2-GSEA", showWarnings = FALSE, recursive = TRUE)
dir.create("3-Visualization", showWarnings = FALSE, recursive = TRUE)

# Clean previous outputs
unlink(c("0-Config", "1-DEG", "2-GSEA", "3-Visualization", "Analysis_summary.txt", "sessionInfo.txt"),
       recursive = TRUE, force = TRUE)

# Re-create after cleaning
dir.create("0-Config", showWarnings = FALSE, recursive = TRUE)
dir.create("1-DEG", showWarnings = FALSE, recursive = TRUE)
dir.create("2-GSEA", showWarnings = FALSE, recursive = TRUE)
dir.create("3-Visualization", showWarnings = FALSE, recursive = TRUE)

custom_gene_sets <- list(
  M1_markers = c("Cd86", "Cd80", "Tnf", "Il1b", "Il6", "Nos2", "Cxcl10", "Cxcl9", "Stat1", "Irf5", "Cd68"),
  M2_markers = c("Mrc1", "Arg1", "Cd163", "Msr1", "Il10", "Tgfb1", "Stat6", "Irf4", "Pparg", "Chil3", "Retnla")
)
KEY_GENES <- c("Tnf", "Il1b", "Il6", "Cxcl10", "Nos2")
RUN_TF_ANALYSIS <- FALSE
RUN_COMPARECLUSTER <- FALSE
validate_batch_design(NULL, DESIGN_FORMULA, sample_names = SAMPLE_NAMES)

writeLines(c(
  "# RNAseq_General demo smoke-test configuration snapshot",
  paste0("# Saved: ", Sys.time()),
  paste0("DEG_PVALUE_COLUMN <- ", deparse(DEG_PVALUE_COLUMN)),
  paste0("DEG_LFC_COLUMN <- ", deparse(DEG_LFC_COLUMN)),
  paste0("GSEA_RANK_COLUMN <- ", deparse(GSEA_RANK_COLUMN)),
  paste0("PAIRWISE_TEST_METHOD <- ", deparse(PAIRWISE_TEST_METHOD)),
  paste0("PAIRWISE_P_ADJUST_METHOD <- ", deparse(PAIRWISE_P_ADJUST_METHOD))
), "0-Config/analysis_config_used.R")

# ---- Preprocessing ----
rawcount <- read_count_table(INPUT_FILE, INPUT_FORMAT)
preprocess_summary_rows <- list()
preprocess_summary_rows[["raw_input"]] <- summarize_raw_input(rawcount, GENE_NAME_COL, BIOTYPE_COL, BIOTYPE_FILTER)

if (!is.null(BIOTYPE_COL) && BIOTYPE_COL %in% colnames(rawcount)) {
  rawcount <- rawcount[rawcount[[BIOTYPE_COL]] == BIOTYPE_FILTER, ]
}
preprocess_summary_rows[["after_biotype_filter"]] <- data.frame(step = "after_biotype_filter_rows", value = nrow(rawcount))

count_col_names <- detect_count_columns(rawcount, GENE_NAME_COL, COUNT_COLS)
validate_sample_design(SAMPLE_NAMES, GROUPS, GROUP_LEVELS, COMPARISONS, count_col_names)

countData <- build_count_matrix(rawcount, GENE_NAME_COL, count_col_names, SAMPLE_NAMES,
                                duplicate_report_file = "./1-DEG/Duplicated_gene_symbols.csv")
validate_count_matrix(countData)
countData_before_gene_filter <- countData

group <- factor(GROUPS, levels = GROUP_LEVELS)
colData <- make_col_data(countData, SAMPLE_NAMES, GROUPS, GROUP_LEVELS)

sample_qc <- calculate_sample_qc(countData, colData, group_col = "condition")
sample_qc <- flag_sample_qc(sample_qc,
                            min_library_size = MIN_LIBRARY_SIZE,
                            min_detected_genes = MIN_DETECTED_GENES,
                            max_zero_fraction = MAX_ZERO_FRACTION,
                            min_median_correlation_z = MIN_MEDIAN_CORRELATION_Z)
write.csv(sample_qc, "./3-Visualization/Sample_QC_metrics.csv", row.names = FALSE)

if (length(SAMPLE_EXCLUDE) > 0) {
  res <- apply_sample_exclusion(countData, colData, SAMPLE_EXCLUDE)
  countData <- res$count_data
  colData <- res$col_data
  group <- colData$condition
}

filter_res <- filter_low_count_genes(countData, group, MIN_COUNT)
countData <- filter_res$count_data
min_replicates <- filter_res$min_replicates

# ---- DESeq2 ----
dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData, colData = colData, design = DESIGN_FORMULA)
vsd <- DESeq2::vst(dds, blind = TRUE)
vsd_mat <- SummarizedExperiment::assay(vsd)

write.csv(data.frame(gene_name = rownames(vsd_mat), vsd_mat, check.names = FALSE),
          "./1-DEG/vsd_matrix.csv", row.names = FALSE)
write.csv(data.frame(sample = rownames(colData), as.data.frame(colData), check.names = FALSE),
          "./1-DEG/colData.csv", row.names = FALSE)

group_colors <- make_group_colors(GROUP_LEVELS)
plot_pca_pdf(vsd, GROUP_LEVELS, group_colors, "./3-Visualization/PCA_plot.pdf")
plot_sample_distance_pdf(vsd, "./3-Visualization/Sample_distance_heatmap.pdf")

dds <- DESeq2::DESeq(dds)
res_list <- extract_deseq2_results(dds, COMPARISONS, condition_col = "condition",
                                   alpha = max(THRESHOLD_GRID$p_cutoff), shrink_type = "ashr")

all_gene_deg_dir <- write_all_gene_deg_results(res_list, outdir = "./1-DEG")
deg_by_threshold <- build_deg_threshold_sets(res_list, THRESHOLD_GRID,
                                             pvalue_column = DEG_PVALUE_COLUMN,
                                             lfc_column = DEG_LFC_COLUMN)
deg_summary <- write_deg_threshold_outputs(deg_by_threshold, outdir = "./1-DEG",
                                           threshold_grid = THRESHOLD_GRID,
                                           pvalue_column = DEG_PVALUE_COLUMN,
                                           lfc_column = DEG_LFC_COLUMN)
write_lfc_strategy_summary(res_list, THRESHOLD_GRID, outdir = "./1-DEG",
                           pvalue_column = DEG_PVALUE_COLUMN)
write_deg_diagnostic_summary(res_list, THRESHOLD_GRID, outdir = "./1-DEG")

# Volcano plot
for (comp_name in names(res_list)) {
  plot_volcano_pdf(res_list[[comp_name]], comp_name = comp_name,
                   pvalue_thresh = 0.05, log2fc_thresh = 1,
                   pvalue_column = DEG_PVALUE_COLUMN, lfc_column = DEG_LFC_COLUMN,
                   filename = paste0("./3-Visualization/Volcano_", comp_name, "_standard.pdf"))
}

# Heatmap. Rasterization is auto-disabled when no cairo device is available, so
# the heatmap renders on minimal/headless systems too (just without raster).
top_genes_mad <- head(order(matrixStats::rowMads(vsd_mat), decreasing = TRUE), 1000)
plot_expression_heatmap_pdf(vsd_mat[top_genes_mad, ],
                            filename = "./3-Visualization/Heatmap_top1000_variable_genes.pdf",
                            title = "Top 1000 Variable Genes",
                            group = as.character(colData[colnames(vsd_mat), "condition"]),
                            group_levels = GROUP_LEVELS, group_colors = group_colors,
                            show_row_names = FALSE, show_column_names = FALSE,
                            cluster_columns = TRUE, width = mm_to_in(183), height = mm_to_in(247))

# ---- Enrichment ----
org_db <- org.Mm.eg.db
org_code <- "mmu"
ranked_lists <- lapply(res_list, ranked_gene_list, rank_column = GSEA_RANK_COLUMN)
background_entrez <- map_symbols_to_entrez(rownames(countData), org_db)
background_universe <- unique(background_entrez$ENTREZID)

ora_threshold <- run_threshold_ora(res_list, THRESHOLD_GRID, org_db, background_universe, org_code,
                                   pvalue_column = DEG_PVALUE_COLUMN, lfc_column = DEG_LFC_COLUMN,
                                   outdir = "./2-GSEA", plotdir = "./3-Visualization")

for (comp_name in names(ranked_lists)) {
  entrezList <- make_entrez_ranked_list(ranked_lists[[comp_name]], org_db)
  if (length(entrezList) > 10) {
    ggo <- run_go_gsea(entrezList, org_db = org_db)
    if (!is.null(ggo) && nrow(as.data.frame(ggo)) > 0) {
      write.csv(as.data.frame(ggo), paste0("./2-GSEA/GSEA_GO_", comp_name, ".csv"), row.names = FALSE)
      plot_gsea_suite_pdf(ggo, paste0("./3-Visualization/GSEA_GO_", comp_name), paste("GSEA GO -", comp_name))
    }
  }
}

# ---- Summary ----
summary_report <- paste0(
  "RNA-seq Analysis Summary (Demo Smoke Test)\n",
  "Species: ", SPECIES, " | Samples: ", ncol(countData), " | Genes: ", nrow(countData), "\n",
  "Comparisons: ", length(COMPARISONS), " | Thresholds: ", nrow(THRESHOLD_GRID), "\n"
)
cat(summary_report)
writeLines(summary_report, "./Analysis_summary.txt")
writeLines(capture.output(sessionInfo()), "./sessionInfo.txt")

# ---- Assertions ----
# Core pipeline files that must always be produced.
expected_files <- c(
  "0-Config/analysis_config_used.R",
  "1-DEG/DEG_threshold_summary.csv",
  "1-DEG/all_genes/DESeq2_all_genes_Treatment_vs_Control.csv",
  "1-DEG/standard/DEG_results_Treatment_vs_Control.csv",
  "1-DEG/DEG_lfc_strategy_summary.csv",
  "1-DEG/DEG_diagnostic_summary.csv",
  "2-GSEA/ORA_threshold_summary.csv",
  "3-Visualization/Sample_distance_heatmap.pdf",
  "3-Visualization/Heatmap_top1000_variable_genes.pdf",
  "Analysis_summary.txt",
  "sessionInfo.txt"
)
# ggplot-based PDFs (PCA, volcano) need a working graphics device. They are
# required on CI (which has cairo) but skipped on minimal/headless systems.
ggplot_pdfs <- c(
  "3-Visualization/PCA_plot.pdf",
  "3-Visualization/Volcano_Treatment_vs_Control_standard.pdf"
)
if (.has_working_cairo()) {
  expected_files <- c(expected_files, ggplot_pdfs)
} else {
  message("cairo device unavailable; not requiring ggplot PDFs: ",
          paste(ggplot_pdfs, collapse = ", "))
}

missing_files <- expected_files[!file.exists(expected_files)]
if (length(missing_files) > 0) {
  stop("Smoke test FAILED. Missing expected output files:\n",
       paste(missing_files, collapse = "\n"))
}

# Enrichment is data-dependent; verify the summary exists and the pipeline attempted ORA.
ora_summary <- utils::read.csv("2-GSEA/ORA_threshold_summary.csv", stringsAsFactors = FALSE)
if (!any(ora_summary$SignificantGenes > 0)) {
  warning("Smoke test: no significant DEGs found in any threshold; enrichment outputs may be empty.")
}

# Sanity check DEG counts are non-negative and total tested genes > 0
stopifnot(nrow(deg_summary) > 0)
stopifnot(all(deg_summary$Total >= 0))
stopifnot(nrow(res_list[[1]]) > 0)

cat("\n========================================\n")
cat("General demo smoke test PASSED.\n")
cat("All expected output files were generated.\n")
cat("========================================\n")

# ---- Topic-template demos ----
# Each demo regenerates its own input data and runs its core pipeline. Run them
# via Rscript so each gets a fresh R session and working directory. The smoke
# test fails if any demo fails.
setwd(repo_root)
demo_scripts <- c(
  "demo_RNAseq_General/run_demo.R",
  "demo_RNAseq_limma_voom/run_demo.R",
  "demo_RNAseq_WGCNA/run_demo.R",
  "demo_RNAseq_TME/regenerate_demo_data.R",
  "demo_RNAseq_TME/run_demo.R",
  "demo_RNAseq_TimeCourse/regenerate_demo_data.R",
  "demo_RNAseq_TimeCourse/run_demo.R",
  "demo_RNAseq_TCGA_GEO/regenerate_demo_data.R",
  "demo_RNAseq_TCGA_GEO/run_demo.R"
)

rscript <- file.path(R.home("bin"), "Rscript")
failed_demos <- character(0)
# Capture each demo's stdout+stderr to a per-demo log so CI artifact upload (or a
# local rerun) shows the exact R error instead of just a non-zero exit status.
demo_log_dir <- file.path(repo_root, "examples", "demo_logs")
dir.create(demo_log_dir, showWarnings = FALSE, recursive = TRUE)
# The demo scripts locate their own directory via --file= and chdir into it, but
# only when Rscript resolves a path containing a separator. Invoke each with its
# repo-root-relative path and force the working directory to the repo root so the
# inherited cwd is predictable regardless of where the General demo left it.
setwd(repo_root)
for (script in demo_scripts) {
  script_rel <- file.path("examples", script)
  script_path <- file.path(repo_root, script_rel)
  if (!file.exists(script_path)) {
    cat("\n>>> SKIP (missing):", script, "\n")
    next
  }
  log_file <- file.path(demo_log_dir, paste0(gsub("/", "__", script), ".log"))
  cat("\n----------------------------------------\n")
  cat(">>> Running:", script, "\n")
  cat("----------------------------------------\n")
  status <- system2(rscript, shQuote(script_rel), stdout = log_file, stderr = log_file)
  # Echo the tail of the demo log so it also appears in the main CI log.
  if (file.exists(log_file)) {
    log_lines <- readLines(log_file, warn = FALSE)
    cat(paste(tail(log_lines, 25), collapse = "\n"), "\n")
  }
  if (!identical(status, 0L)) {
    failed_demos <- c(failed_demos, script)
    cat(">>> FAILED:", script, "(exit status", status, ") — see", log_file, "\n")
  } else {
    cat(">>> OK:", script, "\n")
  }
}

cat("\n========================================\n")
if (length(failed_demos) > 0) {
  cat("Smoke test FAILED. Topic demos that failed:\n")
  cat(paste0("  - ", failed_demos, collapse = "\n"), "\n")
  cat("========================================\n")
  stop("One or more topic-template demos failed.")
}
leaked_rplots <- list.files(repo_root, pattern = "^Rplots\\.pdf$", recursive = TRUE, full.names = TRUE)
if (length(leaked_rplots) > 0) {
  stop("Smoke test FAILED. Plot calls escaped their managed devices and created:\n",
       paste(leaked_rplots, collapse = "\n"))
}
cat("All demo smoke tests PASSED (General + limma-voom + WGCNA + TME + TimeCourse + TCGA-GEO).\n")
cat("========================================\n")
