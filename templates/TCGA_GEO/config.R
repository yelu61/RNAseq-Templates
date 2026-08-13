# =============================================================================
# RNAseq_TCGA_GEO — Analysis Configuration
# =============================================================================
# This is the ONLY file you need to edit for a new project.
# `run_analysis.R` sources this file, then runs the full TCGA/GEO public-data
# pipeline (the same steps as notebooks/RNAseq_TCGA_GEO_Template.ipynb)
# non-interactively: 3-way data acquisition, clinical cleanup, DESeq2
# Tumor-vs-Normal, single-gene expression + survival (KM/Cox), ORA + GSEA.
#
# Usage:
#   Rscript run_analysis.R                  # uses ./config.R
#   Rscript run_analysis.R other_config.R   # uses a different config
#
# Outputs use the numbered run layout:
#   0-Config/ 1-DEG/ 2-GSEA/ 3-Visualization/ 6-Survival/
# (OUTDIR below is the run root; keep "." to write the numbered dirs at the
# project root, matching templates/General/ and templates/Limma_Voom/.)
# =============================================================================

options(stringsAsFactors = FALSE)

# ---- 1.1 Species -------------------------------------------------------------
SPECIES        <- "human"          # "human" or "mouse" (drives OrgDb + KEGG organism)

# ---- 1.2 Data source mode ------------------------------------------------------
# Exactly one acquisition path runs:
#   DOWNLOAD_FROM_GDC = TRUE   -> download via TCGAbiolinks (needs network)
#   else DOWNLOAD_FROM_GEO = TRUE -> download a GEO SeriesMatrix (needs network)
#   else                          -> read the LOCAL_* files below (fully offline)
DOWNLOAD_FROM_GDC <- FALSE
DOWNLOAD_FROM_GEO <- FALSE
GEO_ACCESSION  <- "GSE12345"       # e.g. "GSE15459"; used only when DOWNLOAD_FROM_GEO = TRUE

# ---- 1.3 TCGA download parameters (DOWNLOAD_FROM_GDC = TRUE only) --------------
TCGA_PROJECT      <- "TCGA-STAD"   # e.g. TCGA-COAD, TCGA-BRCA
TCGA_DATA_CATEGORY <- "Transcriptome Profiling"
TCGA_DATA_TYPE    <- "Gene Expression Quantification"
TCGA_WORKFLOW     <- "STAR - Counts"
GDC_COUNTS_ASSAY  <- NULL          # NULL = auto-detect; else assay name (e.g. "unstranded")
GDC_TPM_ASSAY     <- NULL          # NULL = auto-detect; else assay name (e.g. "tpm_unstrand")

# Gene symbol mapping (TCGA download only, where row names are ENSEMBL IDs).
GENE_ID_MAP_FILE  <- NULL          # NULL = derive from SummarizedExperiment rowData;
                                   # else path to file with gene_id, gene_name[, gene_type]

# ---- 1.4 Local file mode (offline; used when both download flags are FALSE) ----
LOCAL_COUNTS_FILE   <- "./0-Data/counts.csv"     # genes x samples; first column = gene symbol
LOCAL_TPM_FILE      <- "./0-Data/tpm.csv"        # optional; if NULL, TPM is not used
LOCAL_CLINICAL_FILE <- "./0-Data/clinical.csv"   # must contain a 'barcode' column
LOCAL_GENE_COLUMN   <- NULL          # NULL = auto-detect; set explicitly for numeric Entrez IDs

# ---- 1.5 Analysis parameters ---------------------------------------------------
MIN_COUNT_PER_SAMPLE_FRAC <- 0.75   # keep genes with count > MIN_COUNT in >= this fraction of samples
MIN_COUNT        <- 1
TUMOR_NORMAL_DESIGN <- TRUE         # run DESeq2 Tumor vs Normal + ORA/GSEA
DEG_LFC_CUTOFF   <- 0.5
DEG_PADJ_CUTOFF  <- 0.05

# ---- 1.6 Survival analysis -------------------------------------------------------
GENES_FOR_SURVIVAL  <- c("MKI67")   # genes tested by median split + quartile split + Cox
CLINICAL_VARS_FOR_KM <- c("ajcc_pathologic_stage")  # clinical vars for KM; numeric vars are median-split
TIME_UNIT        <- "month"         # "day", "month", or "year"

# ---- 1.7 Output root -------------------------------------------------------------
# "." = write the numbered run layout at the project root. Set to a path to nest
# the run under that directory instead.
OUTDIR           <- "."

# ---- 1.8 HTML report ---------------------------------------------------------------
# The shared reports/analysis_report.qmd is oriented to the General (DESeq2)
# layout; wiring is present but off by default so a run does not emit a sparse report.
GENERATE_HTML_REPORT <- FALSE
REPORT_TITLE     <- "TCGA/GEO Tumor-vs-Normal Analysis Report"

# =============================================================================
# Derived values & validation (normally no need to edit below this line)
# =============================================================================
if (!SPECIES %in% c("human", "mouse")) stop("SPECIES must be 'human' or 'mouse'.")
if (isTRUE(DOWNLOAD_FROM_GDC) && isTRUE(DOWNLOAD_FROM_GEO)) {
  stop("DOWNLOAD_FROM_GDC and DOWNLOAD_FROM_GEO cannot both be TRUE; pick one data source.")
}
if (!isTRUE(DOWNLOAD_FROM_GDC) && !isTRUE(DOWNLOAD_FROM_GEO)) {
  if (is.null(LOCAL_COUNTS_FILE) || !nzchar(LOCAL_COUNTS_FILE)) {
    stop("Local-file mode requires LOCAL_COUNTS_FILE.")
  }
  if (is.null(LOCAL_CLINICAL_FILE) || !nzchar(LOCAL_CLINICAL_FILE)) {
    stop("Local-file mode requires LOCAL_CLINICAL_FILE.")
  }
}
if (!TIME_UNIT %in% c("day", "month", "year")) stop("TIME_UNIT must be 'day', 'month', or 'year'.")
if (DEG_PADJ_CUTOFF <= 0 || DEG_PADJ_CUTOFF >= 1) stop("DEG_PADJ_CUTOFF must be in (0, 1).")
if (MIN_COUNT_PER_SAMPLE_FRAC <= 0 || MIN_COUNT_PER_SAMPLE_FRAC > 1) {
  stop("MIN_COUNT_PER_SAMPLE_FRAC must be in (0, 1].")
}

.mode <- if (isTRUE(DOWNLOAD_FROM_GDC)) "GDC download" else if (isTRUE(DOWNLOAD_FROM_GEO)) "GEO download" else "local files"
cat("Configuration loaded: source =", .mode, ", species =", SPECIES,
    ", survival genes =", length(GENES_FOR_SURVIVAL), "\n")
