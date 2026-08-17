# =============================================================================
# RNAseq_General — Analysis Configuration
# =============================================================================
# This is the ONLY file you need to edit for a new project.
# `run_analysis.R` sources this file, then runs the full standard pipeline
# (the same steps as notebooks/RNAseq_General.ipynb) non-interactively.
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
  stringsAsFactors = FALSE
)
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
TME_IOBR_METHODS <- c("estimate", "cibersort", "epic", "xcell")
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

default_threshold_config <- THRESHOLD_GRID[THRESHOLD_GRID$name == DEFAULT_THRESHOLD, ]
DEG_P_CUTOFF  <- default_threshold_config$p_cutoff
DEG_LFC_CUTOFF <- default_threshold_config$log2fc
PADJ_THRESH   <- DEG_P_CUTOFF     # backward-compatible aliases used by TF section
LOG2FC_THRESH <- DEG_LFC_CUTOFF

cat("Configuration loaded:", length(SAMPLE_NAMES), "samples,",
    length(COMPARISONS), "comparisons, default threshold =", DEFAULT_THRESHOLD, "\n")
