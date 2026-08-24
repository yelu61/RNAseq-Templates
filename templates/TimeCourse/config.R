# =============================================================================
# RNAseq_TimeCourse — Analysis Configuration
# =============================================================================
# This is the ONLY file you need to edit for a new project.
# `run_analysis.R` sources this file, then runs the full time-course pipeline
# using the same analysis intent and helper library as the TimeCourse notebook:
# time-point mean aggregation -> Mfuzz soft clustering (+ per-cluster ORA) ->
# time-point-vs-baseline DESeq2.
#
# Usage:
#   Rscript run_analysis.R                  # uses ./config.R
#   Rscript run_analysis.R other_config.R   # uses a different config
#
# Outputs use the numbered run layout: 0-Config/ 1-DEG_Timepoint/
# 3-Visualization/ 5-TimeCourse/ (OUTDIR below is the run root; keep "." to write
# the numbered dirs at the project root, matching templates/General/.)
# =============================================================================

options(stringsAsFactors = FALSE)

# ---- 1.1 Species -------------------------------------------------------------
SPECIES        <- "human"          # "human" or "mouse" (drives the OrgDb for ORA)

# ---- 1.2 Normalized expression input (for Mfuzz clustering) ------------------
EXPR_FILE      <- "./1-DEG/vsd_matrix.csv"   # genes x samples; exported by RNAseq_General
META_FILE      <- "./1-DEG/colData.csv"      # must contain sample, condition, and time columns
GENE_COLUMN    <- NULL             # NULL = use rownames/first column
SAMPLE_COLUMN  <- "sample"
TIME_COLUMN    <- "time"           # numeric or factor time column
GROUP_COLUMN   <- "condition"      # optional treatment/group column
TIME_LEVELS    <- NULL             # e.g. c("Day0", "Day7", "Day14", "Day21"); NULL = infer

# ---- 1.3 Mfuzz time-course soft clustering ------------------------------------
RUN_MFUZZ      <- TRUE
MFUZZ_N_CLUSTERS <- 5
MFUZZ_MIN_ACORE  <- 0.7            # membership threshold for core genes
MFUZZ_SEED       <- 2025           # set.seed for reproducible clustering

# ---- 1.4 Raw count input (for time-point vs baseline DEG) ---------------------
RAW_COUNTS_FILE  <- "./0-Data/raw_counts.tsv"   # genes x samples, raw integer counts
COUNT_META_FILE  <- "./0-Data/metadata.csv"     # must contain sample, time, condition; optional subject
COUNT_GENE_COL   <- "gene_name"
COUNT_SAMPLE_COL <- "sample"
COUNT_BIOTYPE_COL    <- NULL           # set NULL to skip biotype filtering
COUNT_BIOTYPE_FILTER <- "protein_coding"

# ---- 1.5 Optional paired design (repeated measures) ---------------------------
SUBJECT_COL    <- NULL             # e.g. "patient_id"; NULL = unpaired

# ---- 1.6 Time-point vs baseline DEG -------------------------------------------
RUN_TIMEPOINT_DEG <- TRUE
BASELINE_TIME  <- NULL             # NULL = earliest TIME_LEVELS; or specify e.g. "Day0"
DEG_PADJ_CUTOFF <- 0.05
DEG_LFC_CUTOFF  <- 0.5
MIN_COUNT       <- 10              # keep genes with >= this many counts (per-group filter)

# ---- 1.7 Output root ----------------------------------------------------------
# "." = write the numbered run layout (1-DEG_Timepoint/ 3-Visualization/
# 5-TimeCourse/) at the project root. Set to a path to nest the run under a dir.
OUTDIR         <- "."

# ---- 1.8 HTML report -----------------------------------------------------------
# The shared reports/analysis_report.qmd is oriented to the General (DESeq2)
# layout; a time-course-specific report is future work. Wiring is present but
# off by default so a run does not emit a sparse report.
GENERATE_HTML_REPORT <- FALSE
REPORT_TITLE     <- "Time-Course (Mfuzz) Analysis Report"

# ---- 1.9 Run lifecycle -------------------------------------------------------
RUN_ROLE        <- "candidate"
PARENT_RUN_ID   <- NA_character_
RUN_CHANGE_NOTE <- ""
RUN_RETENTION   <- "full"

# =============================================================================
# Derived values & validation (normally no need to edit below this line)
# =============================================================================
if (!SPECIES %in% c("human", "mouse")) stop("SPECIES must be 'human' or 'mouse'.")
if (!is.null(TIME_LEVELS) && length(TIME_LEVELS) < 2) {
  stop("TIME_LEVELS must contain at least 2 time points (or NULL to infer).")
}
if (!is.null(BASELINE_TIME) && !is.null(TIME_LEVELS) && !BASELINE_TIME %in% TIME_LEVELS) {
  stop("BASELINE_TIME '", BASELINE_TIME, "' is not present in TIME_LEVELS.")
}
if (isTRUE(RUN_MFUZZ)) {
  if (MFUZZ_N_CLUSTERS < 2) stop("MFUZZ_N_CLUSTERS must be >= 2.")
  if (MFUZZ_MIN_ACORE <= 0 || MFUZZ_MIN_ACORE > 1) stop("MFUZZ_MIN_ACORE must be in (0, 1].")
}
if (DEG_PADJ_CUTOFF <= 0 || DEG_PADJ_CUTOFF >= 1) stop("DEG_PADJ_CUTOFF must be in (0, 1).")
if (isTRUE(RUN_TIMEPOINT_DEG) && !file.exists(RAW_COUNTS_FILE)) {
  stop("RUN_TIMEPOINT_DEG is TRUE but RAW_COUNTS_FILE was not found: ", RAW_COUNTS_FILE,
       "\nProvide raw integer counts, or set RUN_TIMEPOINT_DEG <- FALSE.")
}

cat("Configuration loaded: time-course analysis, species =", SPECIES,
    "| Mfuzz =", RUN_MFUZZ, "| timepoint DEG =", RUN_TIMEPOINT_DEG, "\n")
