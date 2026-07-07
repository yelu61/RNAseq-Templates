#!/usr/bin/env Rscript
# Validate the WGCNA demo notebook core pipeline.
# Run from repository root: Rscript examples/demo_RNAseq_WGCNA/run_demo.R

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_file) == 0 || this_file == "") this_file <- "examples/demo_RNAseq_WGCNA/run_demo.R"
setwd(dirname(this_file))

options(stringsAsFactors = FALSE)

EXPR_FILE <- "./vsd_matrix.csv"
TRAIT_FILE <- "./colData.csv"
GENE_COLUMN <- NULL
SAMPLE_COLUMN <- "sample"
GROUP_COLUMN <- "condition"
MIN_MAD_QUANTILE <- 0.5
NETWORK_TYPE <- "signed"
POWER_VECTOR <- c(1:10, seq(12, 30, 2))
MIN_MODULE_SIZE <- 30
MERGE_CUT_HEIGHT <- 0.25
TARGET_MODULES <- NULL
OUTDIR <- "RNAseq_WGCNA_Output"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
  library(pheatmap)
})
allowWGCNAThreads()

LIB_DIR <- if (dir.exists("RNAseq_lib")) "RNAseq_lib" else "../../RNAseq_lib"
source(file.path(LIB_DIR, "plot_utils.R"))
theme_set(theme_publication())

expr_raw <- read.csv(EXPR_FILE, check.names = FALSE, colClasses = c("character", rep("numeric", 6)))
if (!is.null(GENE_COLUMN) && GENE_COLUMN %in% colnames(expr_raw)) {
  genes <- expr_raw[[GENE_COLUMN]]
  expr <- as.matrix(expr_raw[, setdiff(colnames(expr_raw), GENE_COLUMN), drop = FALSE])
  rownames(expr) <- genes
} else if (!is.numeric(expr_raw[[1]])) {
  genes <- expr_raw[[1]]
  expr <- as.matrix(expr_raw[, -1, drop = FALSE])
  rownames(expr) <- genes
} else {
  genes <- rownames(expr_raw)
  expr <- as.matrix(expr_raw)
}
mode(expr) <- "numeric"
expr <- expr[!duplicated(rownames(expr)) & rownames(expr) != "", , drop = FALSE]

traits <- read.csv(TRAIT_FILE, check.names = FALSE)
common_samples <- intersect(colnames(expr), traits[[SAMPLE_COLUMN]])
expr <- expr[, common_samples, drop = FALSE]
traits <- traits[match(common_samples, traits[[SAMPLE_COLUMN]]), ]
rownames(traits) <- traits[[SAMPLE_COLUMN]]
stopifnot(all(colnames(expr) == rownames(traits)))

gene_mad <- apply(expr, 1, mad, na.rm = TRUE)
expr <- expr[gene_mad >= quantile(gene_mad, MIN_MAD_QUANTILE, na.rm = TRUE), , drop = FALSE]
datExpr <- t(as.matrix(expr))

gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  traits <- traits[rownames(datExpr), , drop = FALSE]
}

sft <- pickSoftThreshold(datExpr, powerVector = POWER_VECTOR, networkType = NETWORK_TYPE, verbose = 0)
soft_power <- sft$powerEstimate
if (is.na(soft_power)) {
  fit_df <- sft$fitIndices
  soft_power <- fit_df$Power[which.max(fit_df$SFT.R.sq)]
}

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

write.csv(data.frame(gene = colnames(datExpr), module = moduleColors), file.path(OUTDIR, "WGCNA_gene_modules.csv"), row.names = FALSE)
saveRDS(list(net = net, moduleColors = moduleColors, MEs = MEs, datExpr = datExpr, traits = traits), file.path(OUTDIR, "WGCNA_network.rds"))

trait_model <- traits %>% dplyr::select(-dplyr::any_of(SAMPLE_COLUMN)) %>% dplyr::select(dplyr::where(~ is.numeric(.) || is.factor(.) || is.character(.)))
trait_model <- trait_model[, colSums(!is.na(trait_model)) > 0, drop = FALSE]
trait_numeric <- model.matrix(~ . - 1, data = trait_model)
moduleTraitCor <- cor(MEs, trait_numeric, use = "p")
moduleTraitP <- corPvalueStudent(moduleTraitCor, nrow(datExpr))
write.csv(moduleTraitCor, file.path(OUTDIR, "Module_trait_correlation.csv"))
write.csv(moduleTraitP, file.path(OUTDIR, "Module_trait_pvalue.csv"))

gene_module <- data.frame(gene = colnames(datExpr), module = moduleColors)
all_modules <- unique(moduleColors[moduleColors != "grey"])
if (!is.null(TARGET_MODULES)) all_modules <- intersect(all_modules, TARGET_MODULES)

hub_list <- list()
for (mod in all_modules) {
  mod_genes <- gene_module$gene[gene_module$module == mod]
  ME <- MEs[[paste0("ME", mod)]]
  kME <- cor(datExpr[, mod_genes, drop = FALSE], ME, use = "p")
  hub <- data.frame(gene = mod_genes, module = mod, kME = as.numeric(kME)) %>% dplyr::arrange(dplyr::desc(abs(kME)))
  hub_list[[mod]] <- hub
  write.csv(hub, file.path(OUTDIR, paste0("Hub_genes_", mod, ".csv")), row.names = FALSE)
}
write.csv(dplyr::bind_rows(hub_list), file.path(OUTDIR, "Hub_genes_all_modules.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))

cat("\n========================================\n")
cat("WGCNA demo PASSED.\n")
cat("Modules found:", paste(all_modules, collapse = ", "), "\n")
cat("Outputs saved to", OUTDIR, "\n")
cat("========================================\n")
