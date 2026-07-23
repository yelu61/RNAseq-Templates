#!/usr/bin/env Rscript
# Validate the TCGA-GEO demo core pipeline in local-file mode (no network).
# Loads a synthetic TCGA-like cohort: raw counts + TPM + clinical, then runs
# Tumor-vs-Normal DESeq2, single-gene expression, KM survival, and ORA/GSEA.
# Run from repository root: Rscript examples/demo_RNAseq_TCGA_GEO/run_demo.R

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_file) == 0 || this_file == "") this_file <- "examples/demo_RNAseq_TCGA_GEO/run_demo.R"
setwd(dirname(this_file))

options(stringsAsFactors = FALSE)

# ---- Parameters (local-file mode) ----
LOCAL_COUNTS_FILE <- "./0-Data/counts.csv"
LOCAL_TPM_FILE <- "./0-Data/tpm.csv"
LOCAL_CLINICAL_FILE <- "./0-Data/clinical.csv"
LOCAL_GENE_COLUMN <- "gene_name"

MIN_COUNT_PER_SAMPLE_FRAC <- 0.5
MIN_COUNT <- 1
TUMOR_NORMAL_DESIGN <- TRUE
DEG_LFC_CUTOFF <- 0.5
DEG_PADJ_CUTOFF <- 0.05

GENES_FOR_SURVIVAL <- c("MKI67")
CLINICAL_VARS_FOR_KM <- c("ajcc_pathologic_stage")
TIME_UNIT <- "month"

OUTDIR <- "RNAseq_TCGA_GEO_Output"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ---- Setup ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(SummarizedExperiment)
  library(DESeq2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(survival)
  library(survminer)
})

LIB_DIR <- if (dir.exists("RNAseq_lib")) "RNAseq_lib" else "../../RNAseq_lib"
source(file.path(LIB_DIR, "plot_utils.R"))
source(file.path(LIB_DIR, "deg_utils.R"))
source(file.path(LIB_DIR, "enrichment_utils.R"))
source(file.path(LIB_DIR, "tcga_utils.R"))
source(file.path(LIB_DIR, "survival_utils.R"))
source(file.path(LIB_DIR, "data_utils.R"))
theme_set(theme_publication())

# ---- Load local counts / TPM / clinical ----
counts_raw <- read_expression_matrix(LOCAL_COUNTS_FILE, gene_column = LOCAL_GENE_COLUMN)
tpm_raw <- read_expression_matrix(LOCAL_TPM_FILE, gene_column = LOCAL_GENE_COLUMN)
validate_expression_contract(tpm_raw, expected = "tpm")
validate_samples_match(colnames(tpm_raw), colnames(counts_raw), context = "TPM vs counts")
clinical_raw <- read_metadata(LOCAL_CLINICAL_FILE, sample_column = "barcode", required_columns = "barcode")

validate_count_matrix(counts_raw)
cat("Counts:", nrow(counts_raw), "genes x", ncol(counts_raw), "samples\n")

# ---- Tumor/Normal grouping ----
colnames(clinical_raw) <- make.names(colnames(clinical_raw), unique = TRUE)
clinical_raw$tissue_type <- infer_tcga_tumor_normal(clinical_raw$barcode)
validate_samples_match(colnames(counts_raw), clinical_raw$barcode, context = "counts vs clinical")

common_samples <- intersect(colnames(counts_raw), clinical_raw$barcode)
counts_raw <- counts_raw[, common_samples, drop = FALSE]
tpm_raw <- tpm_raw[, common_samples, drop = FALSE]
clinical <- clinical_raw[match(common_samples, clinical_raw$barcode), ]
clinical$condition <- factor(clinical$tissue_type, levels = c("Normal", "Tumor"))
stopifnot(all(c("Normal", "Tumor") %in% clinical$condition))
write.csv(clinical, file.path(OUTDIR, "clinical_clean.csv"), row.names = FALSE)
print(table(clinical$condition, useNA = "ifany"))

# ---- Tumor vs Normal DESeq2 ----
keep_genes <- rowSums(counts_raw > MIN_COUNT) >= ceiling(MIN_COUNT_PER_SAMPLE_FRAC * ncol(counts_raw))
counts_filt <- as.matrix(counts_raw[keep_genes, ])
mode(counts_filt) <- "numeric"

colData <- data.frame(row.names = colnames(counts_filt), condition = clinical$condition)
dds <- DESeqDataSetFromMatrix(countData = counts_filt, colData = colData, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Tumor", "Normal"))
res_df <- as.data.frame(res) |>
  rownames_to_column("gene_name") |>
  arrange(padj)
write.csv(res_df, file.path(OUTDIR, "DESeq2_Tumor_vs_Normal.csv"), row.names = FALSE)

sig_genes <- res_df |>
  filter(!is.na(padj), padj < DEG_PADJ_CUTOFF, abs(log2FoldChange) > DEG_LFC_CUTOFF) |>
  pull(gene_name)
cat("Significant DEGs:", length(sig_genes), "\n")

plot_volcano_pdf(
  res_df, comp_name = "Tumor_vs_Normal",
  pvalue_thresh = DEG_PADJ_CUTOFF, log2fc_thresh = DEG_LFC_CUTOFF,
  filename = file.path(OUTDIR, "Volcano_Tumor_vs_Normal.pdf"),
  pvalue_column = "padj", lfc_column = "log2FoldChange"
)

# ---- Single-gene expression + survival ----
expr_log <- log2(as.matrix(tpm_raw) + 1)
surv_df <- prepare_tcga_survival(clinical)
cat("Survival samples:", nrow(surv_df), " | Events:", sum(surv_df$status), "\n")

for (gene in GENES_FOR_SURVIVAL) {
  if (!gene %in% rownames(expr_log)) { warning("Gene not found: ", gene); next }
  plot_tcga_gene_boxplot_pdf(
    expr_log, gene, clinical$condition,
    filename = file.path(OUTDIR, paste0("Expression_boxplot_", gene, ".pdf")),
    title = paste(gene, "expression")
  )
  surv_df[[gene]] <- as.numeric(expr_log[gene, surv_df$barcode])
  plot_km_by_median_pdf(
    surv_df, value_col = gene,
    filename = file.path(OUTDIR, paste0("KM_", gene, ".pdf")),
    title = paste(gene, "survival"), time_unit = TIME_UNIT
  )
}

# Clinical-variable KM
if (length(CLINICAL_VARS_FOR_KM) > 0) {
  run_clinical_km(surv_df, clinical_df = clinical, var_cols = CLINICAL_VARS_FOR_KM,
                  outdir = OUTDIR, time_unit = TIME_UNIT)
}

# Univariate Cox
surv_vars <- intersect(GENES_FOR_SURVIVAL, colnames(surv_df))
if (length(surv_vars) > 0) {
  uni_cox <- run_univariate_cox(surv_df, vars = surv_vars)
  if (!is.null(uni_cox)) {
    write.csv(uni_cox, file.path(OUTDIR, "univariate_Cox.csv"), row.names = FALSE)
    plot_cox_forest_pdf(uni_cox, filename = file.path(OUTDIR, "univariate_Cox_forest.pdf"),
                        title = "Univariate Cox Regression")
  }
}

# ---- ORA / GSEA ----
deg_list <- genes_for_enrichment(res_df, pvalue_thresh = DEG_PADJ_CUTOFF, log2fc_thresh = DEG_LFC_CUTOFF,
                                 pvalue_column = "padj", lfc_column = "log2FoldChange")
universe <- map_symbols_to_entrez(rownames(counts_raw), org.Hs.eg.db)
ego <- run_go_ora(deg_list$sig, org_db = org.Hs.eg.db, universe = universe$ENTREZID)
if (!is.null(ego)) {
  write.csv(as.data.frame(ego), file.path(OUTDIR, "GO_ORA_Tumor_vs_Normal.csv"), row.names = FALSE)
  plot_enrich_suite_pdf(ego, file.path(OUTDIR, "GO_ORA_Tumor_vs_Normal"), "GO ORA")
}

ranked <- ranked_gene_list(res_df, rank_column = "stat")
entrez_ranked <- make_entrez_ranked_list(ranked, org.Hs.eg.db)
gsea_go <- run_go_gsea(entrez_ranked, org_db = org.Hs.eg.db)
if (!is.null(gsea_go)) {
  write.csv(as.data.frame(gsea_go), file.path(OUTDIR, "GO_GSEA_Tumor_vs_Normal.csv"), row.names = FALSE)
  plot_gsea_suite_pdf(gsea_go, file.path(OUTDIR, "GO_GSEA_Tumor_vs_Normal"), "GO GSEA")
}

save.image(file = file.path(OUTDIR, "TCGA_GEO_workspace.Rdata"))
writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))

# ---- Assertions ----
# Core analysis tables must always be produced.
expected_files <- c(
  file.path(OUTDIR, "clinical_clean.csv"),
  file.path(OUTDIR, "DESeq2_Tumor_vs_Normal.csv"),
  file.path(OUTDIR, "KM_clinical_ajcc_pathologic_stage.pdf"),
  file.path(OUTDIR, "clinical_KM_summary.csv"),
  file.path(OUTDIR, "sessionInfo.txt")
)
missing_files <- expected_files[!file.exists(expected_files)]
if (length(missing_files) > 0) {
  stop("TCGA-GEO demo FAILED. Missing expected outputs:\n", paste(missing_files, collapse = "\n"))
}
# ggplot/EnhancedVolcano-based PDFs need a working cairo device; on minimal
# systems they may be skipped. Only require them when cairo is available.
if (isTRUE(capabilities("cairo"))) {
  pdf_files <- c(
    file.path(OUTDIR, "Volcano_Tumor_vs_Normal.pdf"),
    file.path(OUTDIR, "Expression_boxplot_MKI67.pdf"),
    file.path(OUTDIR, "KM_MKI67.pdf")
  )
  missing_pdf <- pdf_files[!file.exists(pdf_files)]
  if (length(missing_pdf) > 0) {
    warning("TCGA-GEO demo: expected PDF figures were not generated:\n",
            paste(missing_pdf, collapse = "\n"))
  }
}
stopifnot(nrow(res_df) > 0, nrow(surv_df) > 0)

cat("\n========================================\n")
cat("TCGA-GEO demo (local-file mode) PASSED.\n")
cat("DEGs (padj<", DEG_PADJ_CUTOFF, ", |LFC|>", DEG_LFC_CUTOFF, "):", length(sig_genes), "\n")
cat("Outputs saved to", OUTDIR, "\n")
cat("========================================\n")
