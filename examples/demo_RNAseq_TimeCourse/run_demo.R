#!/usr/bin/env Rscript
# Drive the TimeCourse production template against the demo data.
# Thin driver: builds a demo config (absolute paths), then calls
# templates/TimeCourse/run_analysis.R and asserts the expected output files exist.
# Run from repository root: Rscript examples/demo_RNAseq_TimeCourse/run_demo.R

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_file) == 0 || this_file == "") this_file <- "examples/demo_RNAseq_TimeCourse/run_demo.R"
setwd(dirname(this_file))
demo_dir <- normalizePath(getwd())

options(stringsAsFactors = FALSE)

repo_root <- tryCatch(rprojroot::find_root(rprojroot::is_git_root),
                      error = function(e) normalizePath(file.path(demo_dir, "..", "..")))
runner  <- file.path(repo_root, "templates", "TimeCourse", "run_analysis.R")
lib_dir <- file.path(repo_root, "RNAseq_lib")
stopifnot(file.exists(runner), dir.exists(lib_dir))
# NOTE: this macOS machine mangles spaces in system()/system2() arguments to
# "~+~", so an absolute runner path (which contains "Mobile Documents") cannot be
# passed to the child R process. Pass a RELATIVE, space-free runner path instead
# (resolved by the child against the inherited cwd = demo_dir), and hand the
# absolute RNAseq_lib path via the inherited environment.
runner_rel <- file.path("..", "..", "templates", "TimeCourse", "run_analysis.R")

# The orchestrator regenerates inputs before this demo; regenerate here too if a
# required input is missing so the demo is runnable standalone.
required_inputs <- file.path(demo_dir, "0-Data",
                             c("vsd_matrix.csv", "colData.csv", "raw_counts.tsv", "metadata.csv"))
if (any(!file.exists(required_inputs))) {
  regen <- file.path(demo_dir, "regenerate_demo_data.R")
  message("Regenerating demo inputs via ", regen)
  s <- system2("Rscript", shQuote(regen))
  if (!identical(s, 0L) && s != 0) stop("regenerate_demo_data.R failed with status ", s)
}
stopifnot(all(file.exists(required_inputs)))

OUTDIR_NAME <- "RNAseq_TimeCourse_Output"
outdir_abs  <- file.path(demo_dir, OUTDIR_NAME)

# Build a demo config with absolute paths so the child process is independent of
# the caller and the smoke test never writes into a template source directory.
q <- function(x) sprintf('"%s"', x)
config_lines <- c(
  "options(stringsAsFactors = FALSE)",
  'SPECIES        <- "human"',
  paste0("EXPR_FILE      <- ", q(file.path(demo_dir, "0-Data", "vsd_matrix.csv"))),
  paste0("META_FILE      <- ", q(file.path(demo_dir, "0-Data", "colData.csv"))),
  'GENE_COLUMN    <- "gene_name"',
  'SAMPLE_COLUMN  <- "sample"',
  'TIME_COLUMN    <- "time"',
  'GROUP_COLUMN   <- "condition"',
  'TIME_LEVELS    <- c("Day0","Day7","Day14","Day21")',
  "RUN_MFUZZ      <- TRUE",
  "MFUZZ_N_CLUSTERS <- 4",
  "MFUZZ_MIN_ACORE  <- 0.7",
  "MFUZZ_SEED       <- 2025",
  paste0("RAW_COUNTS_FILE  <- ", q(file.path(demo_dir, "0-Data", "raw_counts.tsv"))),
  paste0("COUNT_META_FILE  <- ", q(file.path(demo_dir, "0-Data", "metadata.csv"))),
  'COUNT_GENE_COL   <- "gene_name"',
  'COUNT_SAMPLE_COL <- "sample"',
  "COUNT_BIOTYPE_COL    <- NULL",
  'COUNT_BIOTYPE_FILTER <- "protein_coding"',
  "SUBJECT_COL    <- NULL",
  "RUN_TIMEPOINT_DEG <- TRUE",
  'BASELINE_TIME  <- "Day0"',
  "DEG_PADJ_CUTOFF <- 0.05",
  "DEG_LFC_CUTOFF  <- 0.5",
  "MIN_COUNT       <- 10",
  paste0("OUTDIR         <- ", q(outdir_abs)),
  "GENERATE_HTML_REPORT <- FALSE",
  'REPORT_TITLE     <- "Time-course demo"'
)
config_path <- tempfile(fileext = ".R")
writeLines(config_lines, config_path)

# Clean the OUTDIR so assertions reflect THIS run.
unlink(outdir_abs, recursive = TRUE, force = TRUE)

Sys.setenv(RNASEQ_LIB_DIR = lib_dir)
status <- system2("Rscript", c(shQuote(runner_rel), shQuote(config_path)))
if (is.null(status) || status != 0L) {
  stop("TimeCourse runner exited with status ", status, ".")
}

expected_files <- file.path(outdir_abs, c(
  "0-Config/analysis_config_used.R",
  "5-TimeCourse/mfuzz_clusters.csv",
  "1-DEG_Timepoint/Timepoint_DEG_summary.csv",
  "Analysis_summary.txt",
  "sessionInfo.txt",
  "run_manifest.csv"
))
missing_files <- expected_files[!file.exists(expected_files)]
if (length(missing_files) > 0) {
  stop("TimeCourse demo FAILED (runner exit status ", status, "). Missing expected outputs:\n",
       paste(missing_files, collapse = "\n"))
}

cat("\n========================================\n")
cat("TimeCourse demo PASSED (drove templates/TimeCourse/run_analysis.R).\n")
cat("Outputs saved to", outdir_abs, "\n")
cat("========================================\n")
