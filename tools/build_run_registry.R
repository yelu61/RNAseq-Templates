#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript tools/build_run_registry.R <analysis/runs>")
}

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", cmd, value = TRUE))
file_arg <- gsub("~\\+~", " ", file_arg)
script_dir <- dirname(normalizePath(file_arg[[1]], mustWork = TRUE))
repo_root <- dirname(script_dir)
source(file.path(repo_root, "RNAseq_lib", "run_utils.R"))

runs_dir <- normalizePath(path.expand(args[[1]]), mustWork = TRUE)
registry <- build_run_registry(runs_dir)
cat("Run registry written:", file.path(runs_dir, "RUN_REGISTRY.csv"), "\n")
cat("Runs:", nrow(registry), " Exact duplicates:", sum(nzchar(registry$duplicate_of)), "\n")
