# =============================================================================
# RNAseq_Limma_Voom — Analysis Configuration
# =============================================================================
# This is the ONLY file you need to edit for a new project.
# `run_analysis.R` sources this file, then runs the full limma-voom pipeline
# using the same analysis intent and helper library as the limma-voom notebook.
#
# Usage:
#   Rscript run_analysis.R                  # uses ./config.R
#   Rscript run_analysis.R other_config.R   # uses a different config
#
# Outputs use the numbered run layout: 0-Config/ 1-DEG/ 2-GSEA/ 3-Visualization/
# (OUTDIR below is the run root; keep "." to write the numbered dirs at the
# project root, matching templates/General/.)
# =============================================================================

options(stringsAsFactors = FALSE)

# ---- 1.1 Species -------------------------------------------------------------
SPECIES        <- "human"          # "human" or "mouse" (drives OrgDb + KEGG organism)

# ---- 1.2 Data Input ----------------------------------------------------------
INPUT_FILE     <- "./0-Data/featureCounts_merged_count.annot.tsv"
INPUT_FORMAT   <- "tsv"            # "tsv", "csv", or "excel"
GENE_NAME_COL  <- "gene_name"
BIOTYPE_COL    <- "gene_biotype"   # set NULL to skip biotype filtering
BIOTYPE_FILTER <- "protein_coding"
COUNT_COLS     <- NULL             # NULL = auto-detect numeric columns

# ---- 1.3 Experimental Design -------------------------------------------------
SAMPLE_NAMES   <- c("Control_1", "Control_2", "Control_3",
                    "Treatment_1", "Treatment_2", "Treatment_3")
GROUPS         <- c(rep("Control", 3), rep("Treatment", 3))
GROUP_LEVELS   <- c("Control", "Treatment")   # control level first

# Optional batch covariate (same length as SAMPLE_NAMES); added to the model.
BATCH_VECTOR   <- NULL             # e.g. c("B1","B1","B2","B2","B1","B2")

# ---- 1.4 Comparisons ----------------------------------------------------------
# Each entry: c("output_name", "numerator", "denominator")
COMPARISONS    <- list(
  c("Treatment_vs_Control", "Treatment", "Control")
)

# ---- 1.5 Thresholds -----------------------------------------------------------
DEG_PADJ_CUTOFF  <- 0.05
DEG_LFC_CUTOFF   <- 0.5
MIN_COUNT        <- 10             # keep genes with >= this many counts ...
MIN_SAMPLE_FRAC  <- 0.5            # ... in at least this fraction of samples

# ---- 1.6 Output root -----------------------------------------------------------
# "." = write the numbered run layout (1-DEG/ 2-GSEA/ 3-Visualization/) at the
# project root. Set to a path to nest the run under that directory instead.
OUTDIR         <- "."

# ---- 1.7 HTML report -----------------------------------------------------------
# The shared reports/analysis_report.qmd is oriented to the General (DESeq2)
# layout; a limma-specific report is future work. Wiring is present but off by
# default so a run does not emit a sparse report.
GENERATE_HTML_REPORT <- FALSE
REPORT_TITLE     <- "limma-voom Differential Expression Report"

# ---- 1.8 Run lifecycle -------------------------------------------------------
RUN_ROLE        <- "candidate"
PARENT_RUN_ID   <- NA_character_
RUN_CHANGE_NOTE <- ""
RUN_RETENTION   <- "full"

# =============================================================================
# Derived values & validation (normally no need to edit below this line)
# =============================================================================
if (!SPECIES %in% c("human", "mouse")) stop("SPECIES must be 'human' or 'mouse'.")
if (length(SAMPLE_NAMES) != length(GROUPS)) stop("SAMPLE_NAMES and GROUPS must have equal length.")
if (!all(GROUPS %in% GROUP_LEVELS)) stop("Every GROUPS entry must appear in GROUP_LEVELS.")
if (!is.null(BATCH_VECTOR) && length(BATCH_VECTOR) != length(SAMPLE_NAMES)) {
  stop("BATCH_VECTOR length must match SAMPLE_NAMES.")
}
for (cmp in COMPARISONS) {
  if (length(cmp) != 3) stop("Each COMPARISONS entry must be c(name, numerator, denominator).")
  if (!cmp[2] %in% GROUP_LEVELS || !cmp[3] %in% GROUP_LEVELS) {
    stop("Comparison '", cmp[1], "' references a level not present in GROUP_LEVELS.")
  }
}
if (DEG_PADJ_CUTOFF <= 0 || DEG_PADJ_CUTOFF >= 1) stop("DEG_PADJ_CUTOFF must be in (0, 1).")

cat("Configuration loaded:", length(SAMPLE_NAMES), "samples,",
    length(COMPARISONS), "comparison(s), species =", SPECIES, "\n")
