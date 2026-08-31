#!/usr/bin/env Rscript
# Install all R/Bioconductor dependencies for RNAseq-Templates.
# Run from the repository root: Rscript install_dependencies.R
#
# On CI (GitHub Actions ubuntu + r-lib/actions setup-r with use-public-rspm),
# set P3M as the repo so Linux gets pre-built binaries instead of slow source
# builds. Locally (macOS) the default CRAN repo already serves binaries.

# Prefer Posit Package Manager binaries on Linux; fall back to CRAN elsewhere.
# P3M serves distro-specific binaries via the __linux__/<distro>/ prefix, and the
# distro MUST match the runner's glibc/system-library ABI. ubuntu-latest is now
# noble (24.04); using the jammy binaries pulled a libMagick++ soname (.so.8) that
# noble's apt packages (libmagick++-6.q16-9t64 -> .so.9) do not provide, which
# broke SpatialExperiment/GSVA/IOBR at load time. Detect the codename at runtime
# and fall back to noble.
.on_ci <- nzchar(Sys.getenv("GITHUB_ACTIONS"))
if (.on_ci || (Sys.info()[["sysname"]] == "Linux" && !nzchar(Sys.getenv("RENV_PATHS_ROOT")))) {
  codename <- "noble"
  os_release <- tryCatch(readLines("/etc/os-release", warn = FALSE), error = function(e) character(0))
  vc_line <- grep("^VERSION_CODENAME=", os_release, value = TRUE)
  if (length(vc_line) > 0) {
    codename <- sub('^VERSION_CODENAME="?([^"]*)"?.*$', "\\1", vc_line[1])
  }
  p3m <- sprintf("https://packagemanager.posit.co/cran/__linux__/%s/latest", codename)
  options(repos = c(P3M = p3m))
  message("Using P3M repo for distro '", codename, "': ", p3m)
} else {
  options(repos = c(CRAN = "https://cloud.r-project.org/"))
}
options(
  timeout = 600,            # allow slow CDN downloads of large binaries
  warn = 1,
  Ncpus = max(1L, parallel::detectCores() - 1L)
)
# Binary installs are only meaningful on macOS (CRAN binaries). P3M serves Linux
# binaries through the same install.packages() source path (the repo URL carries
# the binary payload), so a plain "source" type still downloads pre-built .tar.gz
# binaries from P3M. Keep "binary" only for macOS where CRAN uses that type.
.pkg_type <- if (Sys.info()[["sysname"]] == "Darwin") "binary" else "source"
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
  "testthat",       # regression test suite
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
  install.packages(missing_cran, type = .pkg_type)
} else {
  message("All CRAN packages already installed.")
}

# ESTIMATE is distributed through R-Forge rather than CRAN.
if (!requireNamespace("estimate", quietly = TRUE)) {
  message("Installing estimate from R-Forge...")
  install.packages("estimate", repos = c(RForge = "https://R-Forge.R-project.org", getOption("repos")), type = .pkg_type)
}

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", type = .pkg_type)
missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc) > 0) {
  message("Installing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
  # BiocManager resolves Bioc repos itself; it uses getOption("repos") for the CRAN
  # mirror, which we already pointed at P3M on Linux for fast binaries.
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE, type = .pkg_type)
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
  # Hard-fail so CI surfaces the missing packages at the install step instead of
  # "passing with warnings" and only erroring later in the smoke test (which made
  # failures hard to attribute). Local users still get the same clear message.
  stop("Failed to install required packages: ", paste(failed, collapse = ", "),
       "\nInstall the missing system libraries (e.g. libmagick++-dev for ",
       "SpatialExperiment/GSVA) and re-run install_dependencies.R.")
} else {
  message("All dependencies installed successfully.")
}

# Optional: the HTML report uses the quarto CLI when available, else rmarkdown.
if (!nzchar(Sys.which("quarto")) && !requireNamespace("rmarkdown", quietly = TRUE)) {
  message("Note: HTML report generation is disabled. Install the quarto CLI ",
          "(https://quarto.org) or the 'rmarkdown' package to enable it.")
}
