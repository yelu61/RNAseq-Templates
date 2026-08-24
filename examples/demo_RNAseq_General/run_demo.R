#!/usr/bin/env Rscript
# Drive the General production runner against the bundled demo data and verify
# the production-only provenance contract.

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_file) == 0 || this_file == "") this_file <- "examples/demo_RNAseq_General/run_demo.R"
setwd(dirname(this_file))
demo_dir <- normalizePath(getwd())

repo_root <- tryCatch(rprojroot::find_root(rprojroot::is_git_root),
                      error = function(e) normalizePath(file.path(demo_dir, "..", "..")))
lib_dir <- file.path(repo_root, "RNAseq_lib")
base_config <- file.path(repo_root, "templates", "General", "config.R")
input_file <- file.path(repo_root, "examples", "demo_data", "demo_counts.tsv")
stopifnot(dir.exists(lib_dir), file.exists(base_config), file.exists(input_file))

outdir_abs <- file.path(demo_dir, "RNAseq_General_Output")
unlink(outdir_abs, recursive = TRUE, force = TRUE)
dir.create(outdir_abs, recursive = TRUE)

q <- function(x) sprintf('"%s"', x)
config_lines <- c(
  paste0("source(", q(base_config), ")"),
  paste0("INPUT_FILE <- ", q(input_file)),
  'SAMPLE_NAMES <- c("Control_1","Control_2","Control_3","Treatment_1","Treatment_2","Treatment_3")',
  'GROUPS <- c(rep("Control",3), rep("Treatment",3))',
  'GROUP_LEVELS <- c("Control","Treatment")',
  'COMPARISONS <- list(c("Treatment_vs_Control","Treatment","Control"))',
  'custom_gene_sets <- list(M1 = c("Cd86","Tnf","Il1b"), M2 = c("Mrc1","Arg1","Il10"))',
  'KEY_GENES <- c("Tnf","Il1b","Il6")',
  "RUN_COMPARECLUSTER <- FALSE",
  "RUN_TF_ANALYSIS <- FALSE",
  "RUN_TME <- FALSE",
  "EXPORT_EXCEL <- FALSE",
  "GENERATE_HTML_REPORT <- FALSE",
  'RUN_ROLE <- "repro_check"',
  'RUN_CHANGE_NOTE <- "General CLI smoke test"',
  'RUN_RETENTION <- "metadata_only"'
)
config_path <- tempfile(fileext = ".R")
writeLines(config_lines, config_path)

# Invoke from the output root: the production runner deliberately treats its
# invocation directory as the immutable run bundle.
runner_rel <- file.path("..", "..", "..", "templates", "General", "run_analysis.R")
setwd(outdir_abs)
Sys.setenv(RNASEQ_LIB_DIR = lib_dir)
status <- system2("Rscript", c(shQuote(runner_rel), shQuote(config_path)))
setwd(demo_dir)
if (is.null(status) || status != 0L) stop("General runner exited with status ", status, ".")

expected_files <- file.path(outdir_abs, c(
  "0-Config/analysis_config_used.R",
  "1-DEG/DEG_threshold_summary.csv",
  "1-DEG/DEG_results.Rdata",
  "3-Visualization/PCA_plot.pdf",
  "Analysis_summary.txt",
  "sessionInfo.txt",
  "run_manifest.csv"
))
missing_files <- expected_files[!file.exists(expected_files)]
if (length(missing_files) > 0) {
  stop("General CLI demo FAILED. Missing expected outputs:\n", paste(missing_files, collapse = "\n"))
}

manifest <- read.csv(file.path(outdir_abs, "run_manifest.csv"), stringsAsFactors = FALSE)
stopifnot(
  nrow(manifest) == 1,
  manifest$status[[1]] == "completed",
  manifest$role[[1]] == "repro_check",
  grepl("^[0-9a-f]{32}$", manifest$input_md5[[1]]),
  grepl("^[0-9a-f]{32}$", manifest$analysis_signature[[1]]),
  grepl("^[0-9a-f]{32}$", manifest$backend_signature[[1]])
)

cat("\n========================================\n")
cat("General CLI demo PASSED (drove templates/General/run_analysis.R).\n")
cat("Outputs saved to", outdir_abs, "\n")
cat("========================================\n")
