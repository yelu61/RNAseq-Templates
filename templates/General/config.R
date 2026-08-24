# =============================================================================
# RNAseq_General — Analysis Configuration
# =============================================================================
# This is the ONLY file you need to edit for a new project.
# `run_analysis.R` sources this file, then runs the full standard pipeline
# using the same analysis intent and helper library as the General notebook.
#
# Usage:
#   Rscript run_analysis.R                  # uses ./config.R
#   Rscript run_analysis.R other_config.R   # uses a different config
#
# Keep the notebook for interactive work (tweaking gene sets for GSVA, custom
# visualizations); use this script + config to reproduce the standard run.
# =============================================================================

options(stringsAsFactors = FALSE)

# ---- 1.1 Species -------------------------------------------------------------
SPECIES        <- "mouse"          # "mouse" or "human"

# ---- 1.2 Data Input ----------------------------------------------------------
INPUT_FILE     <- "./0-Data/featureCounts_merged_count.annot.tsv"
INPUT_FORMAT   <- "tsv"            # "tsv" or "excel"
GENE_NAME_COL  <- "gene_name"      # gene name column
BIOTYPE_COL    <- "gene_biotype"   # biotype column (set NULL to skip filtering)
BIOTYPE_FILTER <- "protein_coding" # biotype to keep
COUNT_COLS     <- NULL             # NULL = auto-detect; or e.g. c(2:10)

# ---- 1.3 Experimental Design -------------------------------------------------
SAMPLE_NAMES   <- c("EM_LPS_1", "EM_LPS_2", "EM_LPS_3",
                    "LPS_1", "LPS_2", "LPS_3",
                    "blank_1", "blank_2", "blank_3")
GROUPS         <- c(rep("EM_LPS", 3), rep("LPS", 3), rep("blank", 3))
GROUP_LEVELS   <- c("blank", "LPS", "EM_LPS")   # control level first

# Optional batch vector (same length as SAMPLE_NAMES) for batch diagnostics.
BATCH_VECTOR   <- NULL             # e.g. c("B1","B1","B2","B2", ...)

# Optional paired / repeated-measures design.
# If PAIR_ID is provided the design becomes ~ pair_id + condition automatically.
PAIR_ID        <- NULL             # e.g. c("P1","P1","P2","P2","P3","P3")

# ---- 1.4 Sample QC and optional exclusion ------------------------------------
SAMPLE_EXCLUDE <- character(0)     # e.g. c("bad_sample_1")
MIN_LIBRARY_SIZE <- NULL           # e.g. 1e6; NULL = do not flag by this metric
MIN_DETECTED_GENES <- NULL         # e.g. 8000
MAX_ZERO_FRACTION <- NULL          # e.g. 0.80
MIN_MEDIAN_CORRELATION_Z <- -3     # flag if median sample-correlation z < this

# ---- 1.5 Comparisons ----------------------------------------------------------
# Each entry: c("output_name", "treatment", "control")
COMPARISONS    <- list(
  c("LPS_vs_Blank", "LPS",    "blank"),
  c("EM_vs_Blank",  "EM_LPS", "blank"),
  c("EM_vs_LPS",    "EM_LPS", "LPS")
)

# ---- 1.6 Thresholds -----------------------------------------------------------
DEG_PVALUE_COLUMN <- "padj"                 # "padj" (recommended) or "pvalue"
DEG_LFC_COLUMN    <- "log2FoldChange_raw"   # DEG calls/volcano/heatmap/ORA
GSEA_RANK_COLUMN  <- "stat"                 # preranked GSEA rank metric
THRESHOLD_GRID <- data.frame(
  name     = c("strict", "standard", "loose"),
  p_cutoff = c(0.01, 0.05, 0.10),
  log2fc   = c(1.5, 1.0, 0.5),
  p_column = c("padj", "padj", "padj"),
  stringsAsFactors = FALSE
)
# Optional per-row `p_column` overrides DEG_PVALUE_COLUMN for that tier only.
# Low-DEG projects can append an EXPLORATORY nominal-p tier so ORA still has a
# gene list to work with, e.g. add a row:
#   name = "exploratory", p_cutoff = 0.05, log2fc = 0.5, p_column = "pvalue"
# Nominal-p tiers are hypothesis-generation only: never set one as
# DEFAULT_THRESHOLD, and treat any ORA result from it as exploratory until it
# is confirmed at an adjusted-P tier or by rank-based (GSEA/GSVA) evidence.
DEFAULT_THRESHOLD <- "standard"             # must be one of THRESHOLD_GRID$name

MIN_COUNT      <- 10                        # keep genes with count >= this in enough samples
DESIGN_FORMULA <- ~ condition               # edit to ~ batch + condition if applicable
PAIRWISE_TEST_METHOD <- "t.test"            # GSVA / single-gene plots only
PAIRWISE_P_ADJUST_METHOD <- "BH"            # use "none" only for exploratory checks

# ---- 1.7 GSVA: custom gene sets ----------------------------------------------
# Edit freely. Matching to expression rownames is case-insensitive.
custom_gene_sets <- list(
  M1_markers            = c("Cd86", "Cd80", "Tnf", "Il1b", "Il6", "Nos2", "Cxcl10", "Cxcl9", "Stat1", "Irf5", "Cd68"),
  M2_markers            = c("Mrc1", "Arg1", "Cd163", "Msr1", "Il10", "Tgfb1", "Stat6", "Irf4", "Pparg", "Chil3", "Retnla"),
  Inflammatory_response = c("Il1b", "Il6", "Tnf", "Il12b", "Cxcl1", "Cxcl2", "Ccl2", "Ptgs2", "Nfkb1", "Nfkbia"),
  Anti_inflammatory     = c("Il10", "Tgfb1", "Il4ra", "Stat6", "Socs1", "Socs3")
)

# ---- 1.8 Key genes for single-gene plots --------------------------------------
KEY_GENES <- c("Tnf", "Il1b", "Il6", "Cxcl10", "Nos2", "Ptgs2", "Mrc1", "Arg1", "Cd163", "Il10")

# ---- 1.9 Advanced analysis switches -------------------------------------------
RUN_TF_ANALYSIS        <- FALSE   # TF activity (requires dorothea + viper)
RUN_COMPARECLUSTER     <- TRUE    # auto-skips when < 3 groups
COMPARECLUSTER_ONTOLOGY <- "BP"   # "BP", "MF", "CC", or "ALL"

# ---- 1.9b Tumor microenvironment (TME) deconvolution --------------------------
# Deconvolution works on TPM (not the VST matrix used elsewhere), so it needs
# per-gene lengths from your annotation. VST/rlog cannot be inverted to TPM and
# is rejected. Mouse symbols are auto-converted to human orthologs via biomaRt
# (requires network). Set RUN_TME <- FALSE to skip the whole block.
RUN_TME          <- FALSE
TME_GENE_LENGTH_COLUMN <- NULL    # column holding gene length; NULL = use start/end below
TME_GENE_LENGTH_UNIT   <- "bp"    # "bp" or "kb" (only when TME_GENE_LENGTH_COLUMN is set)
TME_GENE_START_COL     <- "gene_start"  # used when TME_GENE_LENGTH_COLUMN is NULL
TME_GENE_END_COL       <- "gene_end"
RUN_TME_ESTIMATE <- TRUE          # native ESTIMATE immune/stromal scores
RUN_TME_IOBR     <- TRUE          # IOBR multi-algorithm deconvolution
# Offline-stable defaults. xcell needs a live biomaRt lookup; opt in only when
# that network dependency is acceptable and its reference retrieval is verified.
TME_IOBR_METHODS <- c("estimate", "cibersort", "epic", "mcpcounter")
TME_IOBR_PERM    <- 1000
RUN_TME_SSGSEA   <- TRUE          # ssGSEA over 28 built-in immune signatures
# Mouse->human ortholog mapping cache (.rds path). The first run queries biomaRt
# online and writes the mapping here; later runs are offline-deterministic.
# NULL = always query online (old behaviour).
TME_ORTHOLOG_CACHE <- NULL        # e.g. file.path(PROJECT_ROOT, "data/processed/mouse_human_orthologs.rds")

# ---- 1.10 Excel export --------------------------------------------------------
EXPORT_EXCEL <- TRUE              # write multi-threshold DEG table to 1-DEG/DEG_results.xlsx

# ---- 1.11 HTML report ----------------------------------------------------------
GENERATE_HTML_REPORT <- TRUE      # render RNAseq_report.html (needs quarto CLI or rmarkdown)
REPORT_TITLE <- "RNA-seq Analysis Report"
SCAFFOLD_REPORT_INTERPRETATION <- FALSE  # first run only: create editable per-section report notes

# ---- 1.12 Run lifecycle ------------------------------------------------------
# These fields describe this run; they do not change the statistical analysis.
# Rebuild analysis/runs/RUN_REGISTRY.csv with tools/build_run_registry.R.
RUN_ROLE        <- "candidate"    # candidate/canonical/sensitivity/repro_check/superseded
PARENT_RUN_ID   <- NA_character_  # source run for a sensitivity or replacement run
RUN_CHANGE_NOTE <- ""             # concise reason this run differs from its parent
RUN_RETENTION   <- "full"         # full/slim/metadata_only; never pruned automatically

# =============================================================================
# Derived values & validation (normally no need to edit below this line)
# =============================================================================
if (!all(c("name", "p_cutoff", "log2fc") %in% colnames(THRESHOLD_GRID))) {
  stop("THRESHOLD_GRID must contain name, p_cutoff, and log2fc columns.")
}
if (!DEG_PVALUE_COLUMN %in% c("padj", "pvalue")) stop("DEG_PVALUE_COLUMN must be 'padj' or 'pvalue'.")
if (!is.character(DEG_LFC_COLUMN) || length(DEG_LFC_COLUMN) != 1) stop("DEG_LFC_COLUMN must be a single column name.")
if (any(duplicated(THRESHOLD_GRID$name))) stop("THRESHOLD_GRID$name must be unique.")
if (!DEFAULT_THRESHOLD %in% THRESHOLD_GRID$name) stop("DEFAULT_THRESHOLD must be one of THRESHOLD_GRID$name.")
if ("p_column" %in% colnames(THRESHOLD_GRID)) {
  bad_col <- !is.na(THRESHOLD_GRID$p_column) & !THRESHOLD_GRID$p_column %in% c("padj", "pvalue")
  if (any(bad_col)) stop("THRESHOLD_GRID$p_column must be 'padj' or 'pvalue' (or NA to inherit DEG_PVALUE_COLUMN).")
  default_row_p_col <- THRESHOLD_GRID$p_column[THRESHOLD_GRID$name == DEFAULT_THRESHOLD]
  if (!is.na(default_row_p_col) && default_row_p_col == "pvalue") {
    warning("DEFAULT_THRESHOLD uses nominal pvalue; nominal-p tiers are exploratory and should not headline the analysis.")
  }
}
if (!RUN_ROLE %in% c("candidate", "canonical", "sensitivity", "repro_check", "superseded")) {
  stop("RUN_ROLE must be candidate/canonical/sensitivity/repro_check/superseded.")
}
if (!RUN_RETENTION %in% c("full", "slim", "metadata_only")) {
  stop("RUN_RETENTION must be full/slim/metadata_only.")
}

default_threshold_config <- THRESHOLD_GRID[THRESHOLD_GRID$name == DEFAULT_THRESHOLD, ]
DEG_P_CUTOFF  <- default_threshold_config$p_cutoff
DEG_LFC_CUTOFF <- default_threshold_config$log2fc
DEFAULT_DEG_PVALUE_COLUMN <- if ("p_column" %in% colnames(default_threshold_config) &&
                                 !is.na(default_threshold_config$p_column[[1]]) &&
                                 nzchar(default_threshold_config$p_column[[1]])) {
  default_threshold_config$p_column[[1]]
} else {
  DEG_PVALUE_COLUMN
}
PADJ_THRESH   <- DEG_P_CUTOFF     # backward-compatible aliases used by TF section
LOG2FC_THRESH <- DEG_LFC_CUTOFF

cat("Configuration loaded:", length(SAMPLE_NAMES), "samples,",
    length(COMPARISONS), "comparisons, default threshold =", DEFAULT_THRESHOLD,
    "(", DEFAULT_DEG_PVALUE_COLUMN, ")\n")
