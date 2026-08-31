#!/usr/bin/env Rscript
# Read-only QC of a saved General run; no analysis, dependencies or outputs are changed.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript audit-gsea-qc.R <General_run_dir>")
run_dir <- normalizePath(args[[1]], mustWork = TRUE)
files <- list.files(file.path(run_dir, "2-GSEA"),
                    pattern = "^GSEA_(GO|KEGG)_.*[.]csv$", full.names = TRUE)
if (!length(files)) stop("No saved GSEA GO/KEGG tables found under: ", run_dir)
cat("Read-only QC of saved General GSEA tables\nRun:", run_dir, "\n\n")
for (file in files) {
  x <- utils::read.csv(file, stringsAsFactors = FALSE)
  needed <- c("ID", "NES", "pvalue", "p.adjust")
  if (!all(needed %in% colnames(x))) stop("Missing GSEA columns in: ", file)
  finite <- is.finite(x$NES) & is.finite(x$pvalue) & is.finite(x$p.adjust)
  identifiable <- !is.na(x$ID) & nzchar(x$ID)
  significant <- finite & identifiable & x$p.adjust < 0.05
  cat("File:", basename(file), "\n")
  print(data.frame(
    rows = nrow(x),
    all_fields_NA = sum(rowSums(!is.na(x)) == 0),
    ID_NA = sum(is.na(x$ID)),
    NES_NA = sum(is.na(x$NES)),
    finite_NES_pvalue_padj = sum(finite),
    identifiable_finite_rows = sum(finite & identifiable),
    identifiable_finite_BH_lt_0_05 = sum(significant),
    min_finite_padj = if (any(is.finite(x$p.adjust))) min(x$p.adjust[is.finite(x$p.adjust)]) else NA_real_
  ), row.names = FALSE)
  cat("Naive base-R p.adjust<0.05 subsetting yields",
      nrow(x[x$p.adjust < 0.05, , drop = FALSE]),
      "rows because NA logical indices create NA rows.\n\n")
}

cat("Versions in the inspecting local environment\n")
cat("R:", R.version.string, "\n")
for (package in c("clusterProfiler", "enrichit", "DOSE")) {
  version <- tryCatch(as.character(utils::packageVersion(package)), error = function(e) "not installed")
  cat(package, version, "\n")
}
cat("\nLocal dependency-source inspection\n")
if (requireNamespace("enrichit", quietly = TRUE)) {
  code <- deparse(enrichit::gsea_gson)
  relevant <- grep("p.adjust <-|pvalueCutoff|gsea_res <- gsea_res\\[", code)
  indices <- sort(unique(unlist(lapply(relevant, function(i) seq.int(max(1L, i - 1L), min(length(code), i + 2L))))))
  cat("Excerpt of installed enrichit::gsea_gson:\n")
  for (i in indices) cat(sprintf("%3d %s\n", i, code[[i]]))
  cat("Observed local implementation subsets using gsea_res$pvalue <= pvalueCutoff without an NA guard.\n")
  cat("This base-R subset turns untestable rows into all-NA rows and filters nominal P, not adjusted P.\n")
}
cat("\nInterpretation limits\n")
cat("These CSVs are nominal-P-filtered exports, not the full tested-pathway inventory.\n")
cat("Do not use exported nrow as the number of significant or testable pathways.\n")
cat("All-NA rows are not nonsignificant pathways; their identifiers are lost in the export.\n")
cat("The finite top-ranked rows still require an explicit prespecified FDR threshold.\n")
cat("No input, result, dependency or repository file was modified by this QC script.\n")
