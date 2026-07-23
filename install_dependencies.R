#!/usr/bin/env Rscript
# Install all R/Bioconductor dependencies for RNAseq-Templates.
# Run from the repository root: Rscript install_dependencies.R

options(
  repos = c(CRAN = "https://cloud.r-project.org/"),
  timeout = 600,            # allow slow CDN downloads of large binaries
  warn = 1,
  Ncpus = max(1L, parallel::detectCores() - 1L)
)
# install_dependencies.R is non-interactive: never prompt for library creation or
# dependency updates, and make BiocManager installs non-interactive too.
options(
  install.packages.check.source = "no",
  BiocManager.ask = FALSE
)
if (interactive() == FALSE) {
  # install.packages() asks to create a missing personal library; suppress that.
  Sys.setenv(R_LIBS_USER = .libPaths()[1L])
}

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
  "ashr",           # DESeq2 lfcShrink type = "ashr"
  "ggrepel",        # volcano plot labels
  "data.table",     # used by notebooks and GEOquery
  "matrixStats",    # rowMads / rowSds in visualization
  "jsonlite",       # notebook JSON validation
  "rprojroot",      # robust repository path discovery
  "rmarkdown",      # HTML report rendering (fallback when quarto CLI is absent)
  "babelgene",      # offline mouse-human ortholog mapping
  "e1071",          # native CIBERSORT support
  "future",         # native CIBERSORT parallel backend
  "furrr"           # native CIBERSORT parallel mapping
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
  "biomaRt",         # mouse-to-human ortholog lookup for TME deconvolution
  "preprocessCore"   # native CIBERSORT quantile normalization
)

# Install CRAN packages. type="binary" prefers pre-built macOS binaries so we
# avoid long source compilation; if a binary is unavailable install.packages()
# falls back to source automatically.
missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) {
  message("Installing CRAN packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran, type = "binary")
} else {
  message("All CRAN packages already installed.")
}

# ESTIMATE is distributed through R-Forge rather than CRAN.
if (!requireNamespace("estimate", quietly = TRUE)) {
  message("Installing estimate from R-Forge...")
  install.packages("estimate", repos = c(RForge = "https://R-Forge.R-project.org", getOption("repos")), type = "binary")
}

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", type = "binary")
missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc) > 0) {
  message("Installing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE, type = "binary")
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
check_pkgs <- c(cran_packages, bioc_packages, "estimate", "IOBR")
failed <- check_pkgs[!vapply(check_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(failed) > 0) {
  warning("Failed to install: ", paste(failed, collapse = ", "))
} else {
  message("All dependencies installed successfully.")
}

# Optional: the HTML report uses the quarto CLI when available, else rmarkdown.
if (!nzchar(Sys.which("quarto")) && !requireNamespace("rmarkdown", quietly = TRUE)) {
  message("Note: HTML report generation is disabled. Install the quarto CLI ",
          "(https://quarto.org) or the 'rmarkdown' package to enable it.")
}
