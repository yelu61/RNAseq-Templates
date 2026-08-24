#!/usr/bin/env Rscript
# Drive the Limma_Voom production template against the demo data.
# Thin driver: builds a demo config (absolute paths), then calls
# templates/Limma_Voom/run_analysis.R and asserts the expected output files exist.
# Run from repository root: Rscript examples/demo_RNAseq_limma_voom/run_demo.R

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_file) == 0 || this_file == "") this_file <- "examples/demo_RNAseq_limma_voom/run_demo.R"
setwd(dirname(this_file))
demo_dir <- normalizePath(getwd())

options(stringsAsFactors = FALSE)

repo_root <- tryCatch(rprojroot::find_root(rprojroot::is_git_root),
                      error = function(e) normalizePath(file.path(demo_dir, "..", "..")))
runner  <- file.path(repo_root, "templates", "Limma_Voom", "run_analysis.R")
lib_dir <- file.path(repo_root, "RNAseq_lib")
stopifnot(file.exists(runner), dir.exists(lib_dir))
# NOTE: this macOS machine mangles spaces in system()/system2() arguments to
# "~+~", so an absolute runner path (which contains "Mobile Documents") cannot be
# passed to the child R process. Pass a RELATIVE, space-free runner path instead
# (resolved by the child against the inherited cwd = demo_dir), and hand the
# absolute RNAseq_lib path via the inherited environment (env vars are not
# subject to the shell-command mangling).
runner_rel <- file.path("..", "..", "templates", "Limma_Voom", "run_analysis.R")

# The demo count table is committed; no regeneration needed.
INPUT_FILE <- file.path(demo_dir, "0-Data", "counts.tsv")
stopifnot(file.exists(INPUT_FILE))

OUTDIR_NAME <- "RNAseq_limma_voom_Output"
outdir_abs  <- file.path(demo_dir, OUTDIR_NAME)

# Build a demo config with absolute paths so the child process is independent of
# the caller and the smoke test never writes into a template source directory.
q <- function(x) sprintf('"%s"', x)
config_lines <- c(
  "options(stringsAsFactors = FALSE)",
  paste0("SPECIES        <- ", q("mouse")),
  paste0("INPUT_FILE     <- ", q(INPUT_FILE)),
  'INPUT_FORMAT   <- "tsv"',
  'GENE_NAME_COL  <- "gene_name"',
  'BIOTYPE_COL    <- "gene_biotype"',
  'BIOTYPE_FILTER <- "protein_coding"',
  "COUNT_COLS     <- NULL",
  'SAMPLE_NAMES   <- c("Control_1","Control_2","Control_3","Treatment_1","Treatment_2","Treatment_3")',
  'GROUPS         <- c(rep("Control",3), rep("Treatment",3))',
  'GROUP_LEVELS   <- c("Control","Treatment")',
  "BATCH_VECTOR   <- NULL",
  'COMPARISONS    <- list(c("Treatment_vs_Control","Treatment","Control"))',
  "DEG_PADJ_CUTOFF  <- 0.05",
  "DEG_LFC_CUTOFF   <- 0.5",
  "MIN_COUNT        <- 10",
  "MIN_SAMPLE_FRAC  <- 0.5",
  paste0("OUTDIR         <- ", q(outdir_abs)),
  "GENERATE_HTML_REPORT <- FALSE",
  'REPORT_TITLE     <- "limma-voom demo"'
)
config_path <- tempfile(fileext = ".R")
writeLines(config_lines, config_path)

# Clean the OUTDIR so assertions reflect THIS run.
unlink(outdir_abs, recursive = TRUE, force = TRUE)

Sys.setenv(RNASEQ_LIB_DIR = lib_dir)
status <- system2("Rscript", c(shQuote(runner_rel), shQuote(config_path)))
if (is.null(status) || status != 0L) {
  stop("limma-voom runner exited with status ", status, ".")
}

expected_files <- file.path(outdir_abs, c(
  "0-Config/analysis_config_used.R",
  "1-DEG/limma_voom_DEG_summary.csv",
  "1-DEG/limma_voom_results.Rdata",
  "3-Visualization/Volcano_Treatment_vs_Control.pdf",
  "Analysis_summary.txt",
  "sessionInfo.txt",
  "run_manifest.csv"
))
missing_files <- expected_files[!file.exists(expected_files)]
if (length(missing_files) > 0) {
  stop("limma-voom demo FAILED (runner exit status ", status, "). Missing expected outputs:\n",
       paste(missing_files, collapse = "\n"))
}

cat("\n========================================\n")
cat("limma-voom demo PASSED (drove templates/Limma_Voom/run_analysis.R).\n")
cat("Outputs saved to", outdir_abs, "\n")
cat("========================================\n")
