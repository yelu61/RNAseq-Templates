# =============================================================================
# RNAseq_TME — Analysis Configuration
# =============================================================================
# This is the ONLY file you need to edit for a new project.
# `run_analysis.R` sources this file, then runs the full TME deconvolution
# pipeline (the same steps as notebooks/RNAseq_TME_Deconvolution_Template.ipynb)
# non-interactively: expression -> ESTIMATE / IOBR / native CIBERSORT -> ssGSEA.
#
# Usage:
#   Rscript run_analysis.R                  # uses ./config.R
#   Rscript run_analysis.R other_config.R   # uses a different config
#
# Outputs use the numbered run layout: 0-Config/ 4-TME/
# (OUTDIR below is the run root; keep "." to write the numbered dirs at the
# project root, matching templates/Limma_Voom/.)
# =============================================================================

options(stringsAsFactors = FALSE)

# ---- 1.1 Data input ----------------------------------------------------------
# Raw integer counts are the default and recommended source; they are converted
# to TPM using gene lengths before deconvolution. "expression" mode reads a
# pre-computed TPM or log2(TPM+1) matrix (VST/rlog is invalid for deconvolution).
INPUT_MODE        <- "raw_counts"          # "raw_counts" (recommended) or "expression"
RAW_COUNTS_FILE   <- "./0-Data/featureCounts_merged_count.annot.tsv"
RAW_COUNTS_FORMAT <- "tsv"                 # "tsv", "csv", or "excel"
EXPR_FILE         <- "./0-Data/TPM_matrix.csv"
EXPR_UNIT         <- "tpm"                 # expression mode: "tpm" or "log2_tpm"

# Metadata file (e.g. exported by the General/Limma_Voom template); must contain
# the sample and group columns named below.
META_FILE         <- "./1-DEG/colData.csv"

# Gene length information for counts-to-TPM conversion (raw_counts mode only).
# Provide either GENE_LENGTH_COLUMN (in bp or kb per GENE_LENGTH_UNIT) or
# GENE_START_COL + GENE_END_COL.
GENE_COLUMN        <- "gene_id"            # stable unique ID preferred for TPM
GENE_LENGTH_COLUMN <- NULL
GENE_LENGTH_UNIT   <- "bp"                 # "bp" or "kb"
GENE_START_COL     <- "gene_start"
GENE_END_COL       <- "gene_end"

# ---- 1.2 Experimental design ---------------------------------------------------
SAMPLE_COLUMN <- "sample"
GROUP_COLUMN  <- "condition"
GROUP_LEVELS  <- NULL                      # NULL = infer from metadata order

# ---- 1.3 Species ---------------------------------------------------------------
# TME methods use human signatures, so mouse data are converted from MGI symbols
# to HGNC symbols via babelgene/biomaRt (needs network). Keep the notebook guard.
SPECIES       <- "human"                   # "human" or "mouse"

# Optional custom group colors. If NULL, make_group_colors(GROUP_LEVELS) is used.
# Example: GROUP_COLORS <- c("Control" = "#6F6F6F", "Treatment" = "#E07B54")
GROUP_COLORS  <- NULL

# ---- 1.4 ESTIMATE (native implementation) ---------------------------------------
# IOBR's estimate method wraps the same estimate package, so native and IOBR
# estimate are essentially identical. Keep TRUE for independent scores.
RUN_ESTIMATE  <- TRUE

# ---- 1.5 IOBR multi-algorithm deconvolution --------------------------------------
# Requires IOBR plus locally cached method reference data for offline runs.
# Keep FALSE for a deterministic offline base run; enable after cache preflight.
RUN_IOBR      <- FALSE
IOBR_METHODS  <- c("estimate", "cibersort", "epic", "xcell")
IOBR_PERM     <- 1000
IOBR_ARRAYS   <- FALSE

# ---- 1.6 Native CIBERSORT (optional, bundled in references/CIBERSORT/) -----------
RUN_CIBERSORT <- FALSE
CIBERSORT_SCRIPT    <- NULL                # NULL = auto references/CIBERSORT/CIBERSORT.R
CIBERSORT_SIGNATURE <- NULL                # NULL = auto LM22 (human) / cibersort_mouse_22.csv (mouse)
CIBERSORT_PERM      <- 1000
CIBERSORT_QN        <- FALSE               # RNA-seq: FALSE; microarray: TRUE

# When both native CIBERSORT and IOBR cibersort run, emit comparison tables and
# PDFs (scatter + Bland-Altman) automatically.
RUN_CIBERSORT_COMPARISON <- TRUE

# ---- 1.7 Output root -------------------------------------------------------------
# "." = write the numbered run layout (0-Config/ 4-TME/) at the project root.
# Set to a path to nest the run under that directory instead.
OUTDIR        <- "."

# ---- 1.8 HTML report -------------------------------------------------------------
# The shared reports/analysis_report.qmd is oriented to the General (DESeq2)
# layout; a TME-specific report is future work. Wiring is present but off by
# default so a run does not emit a sparse report.
GENERATE_HTML_REPORT <- FALSE
REPORT_TITLE         <- "TME Deconvolution Report"

# =============================================================================
# Derived values & validation (normally no need to edit below this line)
# =============================================================================
if (!INPUT_MODE %in% c("raw_counts", "expression")) {
  stop("INPUT_MODE must be 'raw_counts' or 'expression'.")
}
if (INPUT_MODE == "expression" && !EXPR_UNIT %in% c("tpm", "log2_tpm")) {
  stop("EXPR_UNIT must be 'tpm' or 'log2_tpm'. VST/rlog is invalid for TME deconvolution.")
}
if (!SPECIES %in% c("human", "mouse")) stop("SPECIES must be 'human' or 'mouse'.")
if (!GENE_LENGTH_UNIT %in% c("bp", "kb")) stop("GENE_LENGTH_UNIT must be 'bp' or 'kb'.")
if (!is.null(GROUP_COLORS) && is.null(names(GROUP_COLORS))) {
  stop("GROUP_COLORS must be a named vector, e.g. c(Control = '#6F6F6F').")
}

cat("Configuration loaded: INPUT_MODE =", INPUT_MODE, ", species =", SPECIES,
    ", ESTIMATE =", RUN_ESTIMATE, ", IOBR =", RUN_IOBR, ", CIBERSORT =", RUN_CIBERSORT, "\n")
