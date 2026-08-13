# =============================================================================
# RNAseq_WGCNA — Analysis Configuration
# =============================================================================
# This is the ONLY file you need to edit for a new project.
# `run_analysis.R` sources this file, then runs the full WGCNA pipeline
# (the same steps as notebooks/RNAseq_WGCNA_Template.ipynb) non-interactively:
# expression + trait loading -> MAD filtering / sample QC -> soft-threshold
# selection -> blockwise module detection -> module-trait correlation -> hub
# gene export.
#
# Usage:
#   Rscript run_analysis.R                  # uses ./config.R
#   Rscript run_analysis.R other_config.R   # uses a different config
#
# Outputs are written under OUTDIR (default "5-WGCNA") with the config snapshot
# routed to 0-Config/, matching the numbered run layout used by the other
# templates.
# =============================================================================

options(stringsAsFactors = FALSE)

# ---- 1.1 Data Input ----------------------------------------------------------
# EXPR_FILE: genes x samples normalized expression matrix (VST / rlog /
#   log2(TPM+1)). Do NOT use raw counts — WGCNA assumes ~continuous, roughly
#   scale-free expression.
# TRAIT_FILE: sample metadata; must contain SAMPLE_COLUMN. Extra numeric /
#   categorical columns are model-matrixed into traits for correlation.
EXPR_FILE     <- "./1-DEG/vsd_matrix.csv"
TRAIT_FILE    <- "./1-DEG/colData.csv"
GENE_COLUMN   <- NULL              # gene-ID column name, or NULL to auto-use the first column
SAMPLE_COLUMN <- "sample"          # column in TRAIT_FILE holding sample IDs
GROUP_COLUMN  <- "condition"       # grouping column (kept as a trait candidate)

# ---- 1.2 Gene filtering ------------------------------------------------------
MIN_MAD_QUANTILE <- 0.5            # keep genes with MAD >= this quantile (top variable genes)

# ---- 1.3 Network construction -------------------------------------------------
NETWORK_TYPE  <- "signed"          # "signed" (recommended), "unsigned", or "signed hybrid"
POWER_VECTOR  <- c(1:10, seq(12, 30, 2))   # candidate soft-threshold powers
SOFT_POWER    <- NULL              # NULL = auto-pick from pickSoftThreshold; set an integer to override
MIN_MODULE_SIZE   <- 30            # minimum genes per module
MERGE_CUT_HEIGHT  <- 0.25          # dendrogram cut height for merging similar modules

# ---- 1.4 Hub-gene export ------------------------------------------------------
# Restrict hub-gene export to these module colors, or NULL for all modules.
TARGET_MODULES <- NULL             # e.g. c("blue", "turquoise")

# ---- 1.5 Output root -----------------------------------------------------------
# The notebook writes flat into OUTDIR. Routing OUTDIR to "5-WGCNA" groups the
# WGCNA tables/figures under the numbered run layout; the config snapshot goes
# to 0-Config/ regardless.
OUTDIR        <- "5-WGCNA"

# ---- 1.6 HTML report -----------------------------------------------------------
# The shared reports/analysis_report.qmd is oriented to the General (DESeq2)
# layout; a WGCNA-specific report is future work. Wiring is present but off by
# default so a run does not emit a sparse report.
GENERATE_HTML_REPORT <- FALSE
REPORT_TITLE    <- "WGCNA Co-expression Network Report"

# =============================================================================
# Derived values & validation (normally no need to edit below this line)
# =============================================================================
if (!NETWORK_TYPE %in% c("signed", "unsigned", "signed hybrid")) {
  stop("NETWORK_TYPE must be 'signed', 'unsigned', or 'signed hybrid'.")
}
if (!is.numeric(POWER_VECTOR) || length(POWER_VECTOR) < 2 || any(POWER_VECTOR < 1)) {
  stop("POWER_VECTOR must be a numeric vector of >= 2 powers, all >= 1.")
}
if (!is.null(SOFT_POWER) && (!is.numeric(SOFT_POWER) || length(SOFT_POWER) != 1 || SOFT_POWER < 1)) {
  stop("SOFT_POWER must be NULL or a single number >= 1.")
}
if (MIN_MAD_QUANTILE < 0 || MIN_MAD_QUANTILE >= 1) {
  stop("MIN_MAD_QUANTILE must be in [0, 1).")
}
if (!is.numeric(MIN_MODULE_SIZE) || MIN_MODULE_SIZE < 1) {
  stop("MIN_MODULE_SIZE must be a positive integer.")
}
if (MERGE_CUT_HEIGHT <= 0 || MERGE_CUT_HEIGHT >= 1) {
  stop("MERGE_CUT_HEIGHT must be in (0, 1).")
}

cat("Configuration loaded: network =", NETWORK_TYPE,
    "| MAD quantile =", MIN_MAD_QUANTILE,
    "| min module size =", MIN_MODULE_SIZE,
    "| merge cut =", MERGE_CUT_HEIGHT,
    "| outdir =", OUTDIR, "\n")
