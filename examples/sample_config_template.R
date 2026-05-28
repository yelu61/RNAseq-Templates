# Example configuration block for RNAseq_General.ipynb.
# Copy values from here into the notebook parameter cell and edit for a project.

SPECIES <- "mouse" # "mouse" or "human"
INPUT_FILE <- "./0-Data/featureCounts_merged_count.annot.tsv"
INPUT_FORMAT <- "tsv" # "tsv", "csv", or "excel"
GENE_NAME_COL <- "gene_name"
BIOTYPE_COL <- "gene_biotype"
BIOTYPE_FILTER <- "protein_coding"
COUNT_COLS <- NULL

SAMPLE_NAMES <- c(
  "Control_1", "Control_2", "Control_3",
  "Treatment_1", "Treatment_2", "Treatment_3"
)
GROUPS <- c(rep("Control", 3), rep("Treatment", 3))
GROUP_LEVELS <- c("Control", "Treatment")

COMPARISONS <- list(
  c("Treatment_vs_Control", "Treatment", "Control")
)

THRESHOLD_GRID <- data.frame(
  name = c("strict", "standard", "loose"),
  padj = c(0.01, 0.05, 0.10),
  log2fc = c(1.5, 1.0, 0.5),
  stringsAsFactors = FALSE
)
DEFAULT_THRESHOLD <- "standard"

MIN_COUNT <- 10
DESIGN_FORMULA <- ~ condition
PAIRWISE_TEST_METHOD <- "t.test"

custom_gene_sets <- list(
  IFN_response = c("Ifng", "Cxcl9", "Cxcl10", "Stat1", "Irf1"),
  Cytotoxicity = c("Gzma", "Gzmb", "Prf1", "Nkg7", "Gnly")
)

KEY_GENES <- c("Ifng", "Cxcl9", "Cxcl10", "Gzmb", "Prf1")

RUN_TF_ANALYSIS <- FALSE
RUN_COMPARECLUSTER <- TRUE
COMPARECLUSTER_ONTOLOGY <- "BP"
