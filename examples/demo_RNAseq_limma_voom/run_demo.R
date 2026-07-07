#!/usr/bin/env Rscript
# Validate the limma-voom demo notebook core pipeline.
# Run from repository root: Rscript examples/demo_RNAseq_limma_voom/run_demo.R

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_file) == 0 || this_file == "") this_file <- "examples/demo_RNAseq_limma_voom/run_demo.R"
setwd(dirname(this_file))

options(stringsAsFactors = FALSE)

INPUT_FILE <- "./0-Data/counts.tsv"
INPUT_FORMAT <- "tsv"
GENE_NAME_COL <- "gene_name"
BIOTYPE_COL <- "gene_biotype"
BIOTYPE_FILTER <- "protein_coding"
COUNT_COLS <- NULL

SAMPLE_NAMES <- c(
  "Control_1", "Control_2", "Control_3",
  "Treatment_1", "Treatment_2", "Treatment_3"
)
GROUPS <- c(rep("Control", 3), rep("Treatment", 3))
GROUP_LEVELS <- c("Control", "Treatment")

COMPARISONS <- list(
  c("Treatment_vs_Control", "Treatment", "Control")
)

BATCH_VECTOR <- NULL
DEG_PADJ_CUTOFF <- 0.05
DEG_LFC_CUTOFF <- 0.5
MIN_COUNT <- 10
MIN_SAMPLE_FRAC <- 0.5
OUTDIR <- "RNAseq_limma_voom_Output"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTDIR, "1-DEG"), showWarnings = FALSE)
dir.create(file.path(OUTDIR, "2-GSEA"), showWarnings = FALSE)
dir.create(file.path(OUTDIR, "3-Visualization"), showWarnings = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(limma)
  library(sva)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

LIB_DIR <- if (dir.exists("RNAseq_lib")) "RNAseq_lib" else "../../RNAseq_lib"
source(file.path(LIB_DIR, "plot_utils.R"))
source(file.path(LIB_DIR, "io_utils.R"))
source(file.path(LIB_DIR, "deg_utils.R"))
source(file.path(LIB_DIR, "enrichment_utils.R"))
source(file.path(LIB_DIR, "limma_voom_utils.R"))
theme_set(theme_publication())

rawcount <- read_count_table(INPUT_FILE, INPUT_FORMAT)
if (!is.null(BIOTYPE_COL) && BIOTYPE_COL %in% colnames(rawcount)) {
  rawcount <- rawcount[rawcount[[BIOTYPE_COL]] == BIOTYPE_FILTER, ]
}
count_col_names <- detect_count_columns(rawcount, GENE_NAME_COL, COUNT_COLS)
validate_sample_design(SAMPLE_NAMES, GROUPS, GROUP_LEVELS, COMPARISONS, count_col_names)

countData <- build_count_matrix(rawcount, GENE_NAME_COL, count_col_names, SAMPLE_NAMES)
group <- factor(GROUPS, levels = GROUP_LEVELS)
design <- make_group_design(group)

dge <- prepare_dge_for_voom(countData, group = group,
                             min_counts_per_sample = MIN_COUNT,
                             min_sample_frac = MIN_SAMPLE_FRAC)
v <- run_voom(dge, design = design, plot_file = file.path(OUTDIR, "3-Visualization", "voom_mean_variance_trend.pdf"))
res_list <- run_limma_contrasts(v, design, COMPARISONS)
deg_summary <- write_limma_results(res_list, outdir = file.path(OUTDIR, "1-DEG"))

group_colors <- make_group_colors(GROUP_LEVELS)
for (comp_name in names(res_list)) {
  plot_volcano_pdf(
    res_list[[comp_name]], comp_name = comp_name,
    pvalue_thresh = DEG_PADJ_CUTOFF, log2fc_thresh = DEG_LFC_CUTOFF,
    filename = file.path(OUTDIR, "3-Visualization", paste0("Volcano_", comp_name, ".pdf")),
    pvalue_column = "padj", lfc_column = "log2FoldChange"
  )
}
plot_deg_summary_pdf(deg_summary, filename = file.path(OUTDIR, "3-Visualization", "DEG_summary.pdf"))

species_org <- org.Mm.eg.db
organism_code <- "mmu"
universe <- map_symbols_to_entrez(rownames(countData), species_org)$ENTREZID

for (comp_name in names(res_list)) {
  genes <- genes_for_enrichment(
    res_list[[comp_name]],
    pvalue_thresh = DEG_PADJ_CUTOFF, log2fc_thresh = DEG_LFC_CUTOFF,
    pvalue_column = "padj", lfc_column = "log2FoldChange"
  )
  ego <- run_go_ora(genes$sig, org_db = species_org, universe = universe)
  if (!is.null(ego)) {
    write.csv(as.data.frame(ego), file.path(OUTDIR, "2-GSEA", paste0("GO_ORA_", comp_name, ".csv")), row.names = FALSE)
  }
}

save.image(file = file.path(OUTDIR, "limma_voom_workspace.Rdata"))
writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))

cat("\n========================================\n")
cat("limma-voom demo PASSED.\n")
cat("Outputs saved to", OUTDIR, "\n")
cat("========================================\n")
