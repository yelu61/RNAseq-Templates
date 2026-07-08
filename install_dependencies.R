#!/usr/bin/env Rscript
# Install all R/Bioconductor dependencies for RNAseq-Templates.
# Run from the repository root: Rscript install_dependencies.R

options(repos = c(CRAN = "https://cloud.r-project.org/"))

 cran_packages <- c(
  "tidyverse",      # data manipulation + ggplot2
  "devtools",       # for GitHub installs
  "remotes",        # alternative for GitHub installs
  "pheatmap",       # heatmaps
  "ggpubr",         # publication plots + stat_compare_means
  "RColorBrewer",   # color palettes
  "corrplot",       # correlation plots
  "openxlsx",       # Excel input/output
  "UpSetR",         # DEG overlap visualization
  "estimate",       # ESTIMATE scores (R-Forge)
  "ashr",           # DESeq2 lfcShrink type = "ashr"
  "ggrepel",        # volcano plot labels
  "data.table",     # used by notebooks and GEOquery
  "matrixStats"     # rowMads / rowSds in visualization
)

bioc_packages <- c(
  "BiocManager",
  "DESeq2",
  "clusterProfiler",
  "enrichplot",
  "EnhancedVolcano",
  "GSVA",
  "limma",
  "edgeR",
  "sva",
  "WGCNA",
  "ComplexHeatmap",
  "circlize",
  "SummarizedExperiment",
  "TCGAbiolinks",
  "Mfuzz",
  "survival",
  "survminer",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "msigdbr",        # MSigDB gene sets
  "DOSE",           # enrichment visualization dependency
  "BiocParallel",   # parallel backends
  "GEOquery",        # GEO SeriesMatrix download (optional but pre-installed)
  "biomaRt"          # mouse-to-human ortholog lookup for TME deconvolution
)

# Install CRAN packages
missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) {
  message("Installing CRAN packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran)
} else {
  message("All CRAN packages already installed.")
}

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc) > 0) {
  message("Installing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, ask = FALSE)
} else {
  message("All Bioconductor packages already installed.")
}

# Install IOBR from GitHub (not on CRAN/Bioconductor)
if (!requireNamespace("IOBR", quietly = TRUE)) {
  message("Installing IOBR from GitHub...")
  remotes::install_github("IOBR/IOBR")
} else {
  message("IOBR already installed.")
}

# Verify
check_pkgs <- c(cran_packages, bioc_packages, "IOBR")
failed <- check_pkgs[!vapply(check_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(failed) > 0) {
  warning("Failed to install: ", paste(failed, collapse = ", "))
} else {
  message("All dependencies installed successfully.")
}
