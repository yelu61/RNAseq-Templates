#!/usr/bin/env Rscript
# Drive the WGCNA production template against the demo data.
# Thin driver: builds a demo config (absolute paths), then calls
# templates/WGCNA/run_analysis.R and asserts the expected output files exist.
# Run from repository root: Rscript examples/demo_RNAseq_WGCNA/run_demo.R

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_file) == 0 || this_file == "") this_file <- "examples/demo_RNAseq_WGCNA/run_demo.R"
setwd(dirname(this_file))
demo_dir <- normalizePath(getwd())

options(stringsAsFactors = FALSE)

repo_root <- tryCatch(rprojroot::find_root(rprojroot::is_git_root),
                      error = function(e) normalizePath(file.path(demo_dir, "..", "..")))
runner  <- file.path(repo_root, "templates", "WGCNA", "run_analysis.R")
lib_dir <- file.path(repo_root, "RNAseq_lib")
stopifnot(file.exists(runner), dir.exists(lib_dir))
# NOTE: this macOS machine mangles spaces in system()/system2() arguments to
# "~+~", so an absolute runner path (which contains "Mobile Documents") cannot be
# passed to the child R process. Pass a RELATIVE, space-free runner path instead
# (resolved by the child against the inherited cwd = demo_dir), and hand the
# absolute RNAseq_lib path via the inherited environment.
runner_rel <- file.path("..", "..", "templates", "WGCNA", "run_analysis.R")

# The demo expression + trait tables are committed; no regeneration needed.
EXPR_FILE  <- file.path(demo_dir, "vsd_matrix.csv")
TRAIT_FILE <- file.path(demo_dir, "colData.csv")
stopifnot(file.exists(EXPR_FILE), file.exists(TRAIT_FILE))

# WGCNA run_analysis.R uses RUN_ROOT <- dirname(OUTDIR): tables/figures land in
# OUTDIR (set to .../RNAseq_WGCNA_Output/5-WGCNA) and the run-root artifacts
# (0-Config/, Analysis_summary.txt, sessionInfo.txt) land in its parent.
OUTDIR_NAME <- "RNAseq_WGCNA_Output"
run_root    <- file.path(demo_dir, OUTDIR_NAME)
outdir_abs  <- file.path(run_root, "5-WGCNA")

# Build a demo config with absolute paths so the child process is independent of
# the caller and the smoke test never writes into a template source directory.
q <- function(x) sprintf('"%s"', x)
config_lines <- c(
  "options(stringsAsFactors = FALSE)",
  paste0("EXPR_FILE     <- ", q(EXPR_FILE)),
  paste0("TRAIT_FILE    <- ", q(TRAIT_FILE)),
  'GENE_COLUMN   <- "gene_name"',
  'SAMPLE_COLUMN <- "sample"',
  'GROUP_COLUMN  <- "condition"',
  "MIN_MAD_QUANTILE <- 0.5",
  'NETWORK_TYPE  <- "signed"',
  "POWER_VECTOR  <- c(1:10, seq(12, 30, 2))",
  "SOFT_POWER    <- NULL",
  "MIN_MODULE_SIZE   <- 30",
  "MERGE_CUT_HEIGHT  <- 0.25",
  "TARGET_MODULES <- NULL",
  paste0("OUTDIR        <- ", q(outdir_abs)),
  "GENERATE_HTML_REPORT <- FALSE",
  'REPORT_TITLE    <- "WGCNA demo"'
)
config_path <- tempfile(fileext = ".R")
writeLines(config_lines, config_path)

# Clean the run root so assertions reflect THIS run.
unlink(run_root, recursive = TRUE, force = TRUE)

Sys.setenv(RNASEQ_LIB_DIR = lib_dir)
status <- system2("Rscript", c(shQuote(runner_rel), shQuote(config_path)))
if (is.null(status) || status != 0L) {
  stop("WGCNA runner exited with status ", status, ".")
}

expected_files <- c(
  file.path(run_root, "0-Config", "analysis_config_used.R"),
  file.path(outdir_abs, "WGCNA_network.rds"),
  file.path(outdir_abs, "WGCNA_gene_modules.csv"),
  file.path(outdir_abs, "Module_trait_correlation.csv"),
  file.path(run_root, "Analysis_summary.txt"),
  file.path(run_root, "sessionInfo.txt"),
  file.path(run_root, "run_manifest.csv")
)
missing_files <- expected_files[!file.exists(expected_files)]
if (length(missing_files) > 0) {
  stop("WGCNA demo FAILED (runner exit status ", status, "). Missing expected outputs:\n",
       paste(missing_files, collapse = "\n"))
}

cat("\n========================================\n")
cat("WGCNA demo PASSED (drove templates/WGCNA/run_analysis.R).\n")
cat("Outputs saved to", run_root, "\n")
cat("========================================\n")
