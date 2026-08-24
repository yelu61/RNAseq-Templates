#!/usr/bin/env Rscript
# Drive the TCGA_GEO production template against the demo data (local-file mode).
# Thin driver: builds a demo config (absolute paths), then calls
# templates/TCGA_GEO/run_analysis.R and asserts the expected output files exist.
# Run from repository root: Rscript examples/demo_RNAseq_TCGA_GEO/run_demo.R

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_file) == 0 || this_file == "") this_file <- "examples/demo_RNAseq_TCGA_GEO/run_demo.R"
setwd(dirname(this_file))
demo_dir <- normalizePath(getwd())

options(stringsAsFactors = FALSE)

repo_root <- tryCatch(rprojroot::find_root(rprojroot::is_git_root),
                      error = function(e) normalizePath(file.path(demo_dir, "..", "..")))
runner  <- file.path(repo_root, "templates", "TCGA_GEO", "run_analysis.R")
lib_dir <- file.path(repo_root, "RNAseq_lib")
stopifnot(file.exists(runner), dir.exists(lib_dir))
# NOTE: this macOS machine mangles spaces in system()/system2() arguments to
# "~+~", so an absolute runner path (which contains "Mobile Documents") cannot be
# passed to the child R process. Pass a RELATIVE, space-free runner path instead
# (resolved by the child against the inherited cwd = demo_dir), and hand the
# absolute RNAseq_lib path via the inherited environment.
runner_rel <- file.path("..", "..", "templates", "TCGA_GEO", "run_analysis.R")

# The orchestrator regenerates inputs before this demo; regenerate here too if a
# required input is missing so the demo is runnable standalone.
required_inputs <- file.path(demo_dir, "0-Data", c("counts.csv", "tpm.csv", "clinical.csv"))
if (any(!file.exists(required_inputs))) {
  regen <- file.path(demo_dir, "regenerate_demo_data.R")
  message("Regenerating demo inputs via ", regen)
  s <- system2("Rscript", shQuote(regen))
  if (!identical(s, 0L) && s != 0) stop("regenerate_demo_data.R failed with status ", s)
}
stopifnot(all(file.exists(required_inputs)))

OUTDIR_NAME <- "RNAseq_TCGA_GEO_Output"
outdir_abs  <- file.path(demo_dir, OUTDIR_NAME)

# Build a demo config with absolute paths so the child process is independent of
# the caller and the smoke test never writes into a template source directory. Local-file mode keeps
# the run fully offline.
q <- function(x) sprintf('"%s"', x)
config_lines <- c(
  "options(stringsAsFactors = FALSE)",
  'SPECIES        <- "human"',
  "DOWNLOAD_FROM_GDC <- FALSE",
  "DOWNLOAD_FROM_GEO <- FALSE",
  'GEO_ACCESSION  <- "GSE12345"',
  paste0("LOCAL_COUNTS_FILE   <- ", q(file.path(demo_dir, "0-Data", "counts.csv"))),
  paste0("LOCAL_TPM_FILE      <- ", q(file.path(demo_dir, "0-Data", "tpm.csv"))),
  paste0("LOCAL_CLINICAL_FILE <- ", q(file.path(demo_dir, "0-Data", "clinical.csv"))),
  'LOCAL_GENE_COLUMN   <- "gene_name"',
  "MIN_COUNT_PER_SAMPLE_FRAC <- 0.5",
  "MIN_COUNT        <- 1",
  "TUMOR_NORMAL_DESIGN <- TRUE",
  "DEG_LFC_CUTOFF   <- 0.5",
  "DEG_PADJ_CUTOFF  <- 0.05",
  'GENES_FOR_SURVIVAL  <- c("MKI67")',
  'CLINICAL_VARS_FOR_KM <- c("ajcc_pathologic_stage")',
  'TIME_UNIT        <- "month"',
  paste0("OUTDIR           <- ", q(outdir_abs)),
  "GENERATE_HTML_REPORT <- FALSE",
  'REPORT_TITLE     <- "TCGA-GEO demo"'
)
config_path <- tempfile(fileext = ".R")
writeLines(config_lines, config_path)

# Clean the OUTDIR so assertions reflect THIS run.
unlink(outdir_abs, recursive = TRUE, force = TRUE)

Sys.setenv(RNASEQ_LIB_DIR = lib_dir)
status <- system2("Rscript", c(shQuote(runner_rel), shQuote(config_path)))
if (is.null(status) || status != 0L) {
  stop("TCGA-GEO runner exited with status ", status, ".")
}

expected_files <- file.path(outdir_abs, c(
  "0-Config/analysis_config_used.R",
  "1-DEG/DESeq2_Tumor_vs_Normal.csv",
  "6-Survival/univariate_Cox.csv",
  "3-Visualization/Volcano_Tumor_vs_Normal.pdf",
  "Analysis_summary.txt",
  "sessionInfo.txt",
  "run_manifest.csv"
))
missing_files <- expected_files[!file.exists(expected_files)]
if (length(missing_files) > 0) {
  stop("TCGA-GEO demo FAILED (runner exit status ", status, "). Missing expected outputs:\n",
       paste(missing_files, collapse = "\n"))
}

cat("\n========================================\n")
cat("TCGA-GEO demo (local-file mode) PASSED (drove templates/TCGA_GEO/run_analysis.R).\n")
cat("Outputs saved to", outdir_abs, "\n")
cat("========================================\n")
