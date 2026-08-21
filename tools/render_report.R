#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  stop("Usage: Rscript tools/render_report.R <run_dir> [--title=...] [--author=...] ",
       "[--primary-threshold=standard] [--scaffold] [--publish]")
}

run_dir <- normalizePath(path.expand(args[[1]]), mustWork = TRUE)
options <- args[-1]
value_of <- function(prefix, default) {
  hit <- options[startsWith(options, prefix)]
  if (length(hit)) substring(hit[[1]], nchar(prefix) + 1L) else default
}
title <- value_of("--title=", "RNA-seq Analysis Report")
author <- value_of("--author=", Sys.info()[["user"]])
primary <- value_of("--primary-threshold=", "standard")
scaffold <- "--scaffold" %in% options
publish <- "--publish" %in% options

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", cmd, value = TRUE))
file_arg <- gsub("~\\+~", " ", file_arg)
script_dir <- dirname(normalizePath(file_arg[[1]], mustWork = TRUE))
repo_root <- dirname(script_dir)
source(file.path(repo_root, "RNAseq_lib", "report_utils.R"))

if (scaffold) {
  scaffold_report_interpretation(run_dir)
  scaffold_report_review(run_dir)
}
report <- render_analysis_report(
  outdir = run_dir,
  report_file = file.path(run_dir, "RNAseq_report.html"),
  template = file.path(repo_root, "reports", "analysis_report.qmd"),
  params = list(title = title, author = author, primary_threshold = primary)
)
validate_analysis_report(run_dir, report, publish = publish)
cat(if (publish) "Publish-ready report: " else "Draft report rendered: ", report, "\n", sep = "")
