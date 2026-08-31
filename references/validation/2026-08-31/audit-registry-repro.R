#!/usr/bin/env Rscript
# Minimal v0.11.0 registry reproductions; creates and removes temporary inputs only.
# Run from the repository root, or pass its absolute path as the first argument.
args <- commandArgs(trailingOnly = TRUE)
repo_root <- normalizePath(if (length(args)) args[[1]] else getwd(), mustWork = TRUE)
source(file.path(repo_root, "RNAseq_lib", "run_utils.R"))

main <- function() {
  base <- tempfile("registry-freeze-audit-")
  dir.create(base)
  on.exit(unlink(base, recursive = TRUE), add = TRUE)
  counts <- file.path(base, "counts.csv")
  meta <- file.path(base, "metadata.csv")
  config <- file.path(base, "config.R")
  writeLines("gene,S1\nA,1", counts)
  writeLines("same configuration", config)

  runs <- file.path(base, "metadata_runs")
  dir.create(runs)
  metadata_md5 <- character(0)
  for (i in 1:2) {
    id <- c("before", "after")[[i]]
    run <- file.path(runs, id)
    dir.create(run)
    writeLines(paste0("sample,condition\nS1,", c("control", "treatment")[[i]]), meta)
    metadata_md5[[id]] <- unname(tools::md5sum(meta))
    write_run_manifest(
      run, config, counts, list(META_FILE = meta),
      completed_at = as.POSIXct(i, origin = "2026-01-01", tz = "UTC")
    )
  }
  registry <- build_run_registry(runs)
  cat("CASE 1: auxiliary metadata changes at the same path\n")
  print(metadata_md5)
  print(registry[, c("run_id", "input_md5", "analysis_signature", "duplicate_of")], row.names = FALSE)
  stopifnot(metadata_md5[["before"]] != metadata_md5[["after"]])
  stopifnot(identical(registry$duplicate_of[registry$run_id == "after"], "before"))
  cat("REPRODUCED: changed scientific metadata is nevertheless marked duplicate.\n\n")

  online_runs <- file.path(base, "online_runs")
  dir.create(online_runs)
  for (i in 1:2) {
    run <- file.path(online_runs, paste0("online", i))
    dir.create(run)
    write_run_manifest(
      run, config, NA_character_, list(GEO_ACCESSION = "same"),
      completed_at = as.POSIXct(i, origin = "2026-01-01", tz = "UTC")
    )
  }
  registry <- build_run_registry(online_runs)
  cat("CASE 2: online runs with absent primary-input checksum\n")
  print(registry[, c("run_id", "input_md5", "duplicate_of")], row.names = FALSE)
  cat("R nzchar(NA):", nzchar(NA_character_), "\n")
  stopifnot(all(is.na(registry$input_md5)))
  stopifnot(identical(registry$duplicate_of[registry$run_id == "online2"], "online1"))
  cat("REPRODUCED: missing input checksums are nevertheless grouped as duplicate.\n\n")
  cat("No repository files were changed. Temporary data will be removed.\n")
  print(sessionInfo())
}

main()
