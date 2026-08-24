#!/usr/bin/env Rscript
# Drive the TME production template against the demo data.
# Thin driver: builds a demo config (absolute paths), then calls
# templates/TME/run_analysis.R and asserts the expected output files exist.
# Run from repository root: Rscript examples/demo_RNAseq_TME/run_demo.R

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_file) == 0 || this_file == "") this_file <- "examples/demo_RNAseq_TME/run_demo.R"
setwd(dirname(this_file))
demo_dir <- normalizePath(getwd())

options(stringsAsFactors = FALSE)

repo_root <- tryCatch(rprojroot::find_root(rprojroot::is_git_root),
                      error = function(e) normalizePath(file.path(demo_dir, "..", "..")))
runner  <- file.path(repo_root, "templates", "TME", "run_analysis.R")
lib_dir <- file.path(repo_root, "RNAseq_lib")
stopifnot(file.exists(runner), dir.exists(lib_dir))
# NOTE: this macOS machine mangles spaces in system()/system2() arguments to
# "~+~", so an absolute runner path (which contains "Mobile Documents") cannot be
# passed to the child R process. Pass a RELATIVE, space-free runner path instead
# (resolved by the child against the inherited cwd = demo_dir), and hand the
# absolute RNAseq_lib path via the inherited environment.
runner_rel <- file.path("..", "..", "templates", "TME", "run_analysis.R")

# The orchestrator regenerates inputs before this demo; regenerate here too if a
# required input is missing so the demo is runnable standalone.
required_inputs <- file.path(demo_dir, "0-Data", c("counts.tsv", "metadata.csv"))
if (any(!file.exists(required_inputs))) {
  regen <- file.path(demo_dir, "regenerate_demo_data.R")
  message("Regenerating demo inputs via ", regen)
  s <- system2("Rscript", shQuote(regen))
  if (!identical(s, 0L) && s != 0) stop("regenerate_demo_data.R failed with status ", s)
}
stopifnot(all(file.exists(required_inputs)))

OUTDIR_NAME <- "RNAseq_TME_Deconvolution_Output"
outdir_abs  <- file.path(demo_dir, OUTDIR_NAME)

# Build a demo config with absolute paths so the child process is independent of
# the caller and the smoke test never writes into a template source directory.
q <- function(x) sprintf('"%s"', x)
config_lines <- c(
  "options(stringsAsFactors = FALSE)",
  'INPUT_MODE        <- "raw_counts"',
  paste0("RAW_COUNTS_FILE   <- ", q(file.path(demo_dir, "0-Data", "counts.tsv"))),
  'RAW_COUNTS_FORMAT <- "tsv"',
  paste0("EXPR_FILE         <- ", q(file.path(demo_dir, "0-Data", "TPM_matrix.csv"))),
  'EXPR_UNIT         <- "tpm"',
  paste0("META_FILE         <- ", q(file.path(demo_dir, "0-Data", "metadata.csv"))),
  'GENE_COLUMN        <- "gene_name"',
  "GENE_LENGTH_COLUMN <- NULL",
  'GENE_LENGTH_UNIT   <- "bp"',
  'GENE_START_COL     <- "gene_start"',
  'GENE_END_COL       <- "gene_end"',
  'SAMPLE_COLUMN <- "sample"',
  'GROUP_COLUMN  <- "condition"',
  'GROUP_LEVELS  <- c("Control","Treatment")',
  'SPECIES       <- "human"',
  "GROUP_COLORS  <- NULL",
  "RUN_ESTIMATE  <- TRUE",
  # Native ESTIMATE + ssGSEA give the smoke test deterministic offline TME
  # coverage. IOBR reference bundles are exercised separately when cached.
  "RUN_IOBR      <- FALSE",
  'IOBR_METHODS  <- c("estimate","epic","xcell","cibersort")',
  "IOBR_PERM     <- 100",
  "IOBR_ARRAYS   <- FALSE",
  # Native CIBERSORT needs the external script and is slow; leave it off.
  "RUN_CIBERSORT <- FALSE",
  "CIBERSORT_SCRIPT    <- NULL",
  "CIBERSORT_SIGNATURE <- NULL",
  "CIBERSORT_PERM      <- 100",
  "CIBERSORT_QN        <- FALSE",
  "RUN_CIBERSORT_COMPARISON <- FALSE",
  paste0("OUTDIR        <- ", q(outdir_abs)),
  "GENERATE_HTML_REPORT <- FALSE",
  'REPORT_TITLE         <- "TME demo"'
)
config_path <- tempfile(fileext = ".R")
writeLines(config_lines, config_path)

# Clean the OUTDIR so assertions reflect THIS run.
unlink(outdir_abs, recursive = TRUE, force = TRUE)

Sys.setenv(RNASEQ_LIB_DIR = lib_dir)
status <- system2("Rscript", c(shQuote(runner_rel), shQuote(config_path)))
if (is.null(status) || status != 0L) {
  stop("TME runner exited with status ", status, ".")
}

expected_files <- file.path(outdir_abs, c(
  "0-Config/analysis_config_used.R",
  "4-TME/TPM_matrix.csv",
  "4-TME/ESTIMATE_scores.csv",
  "4-TME/ssGSEA_immune_scores.csv",
  "4-TME/tme_results.Rdata",
  "Analysis_summary.txt",
  "sessionInfo.txt",
  "run_manifest.csv"
))
missing_files <- expected_files[!file.exists(expected_files)]
if (length(missing_files) > 0) {
  stop("TME demo FAILED (runner exit status ", status, "). Missing expected outputs:\n",
       paste(missing_files, collapse = "\n"))
}

cat("\n========================================\n")
cat("TME deconvolution demo PASSED (drove templates/TME/run_analysis.R).\n")
cat("Outputs saved to", outdir_abs, "\n")
cat("========================================\n")
