# Minimal reproduction of v0.11.0 TME mouse/native-CIBERSORT ID mismatch.
# Usage: Rscript audit-tme-repro.R /absolute/path/to/RNAseq-Templates
# Uses installed babelgene data only; no downloads or deconvolution are run.
args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(if (length(args)) args[[1]] else getwd(), mustWork = TRUE)
source(file.path(repo, "RNAseq_lib", "data_utils.R"))
source(file.path(repo, "RNAseq_lib", "tme_utils.R"))
cat(R.version.string, "\nbabelgene", as.character(packageVersion("babelgene")), "\n")
sig <- read_cibersort_signature(file.path(repo, "references", "CIBERSORT", "cibersort_mouse_22.csv"))
expr <- as.data.frame(cbind(S1 = seq_len(nrow(sig)), S2 = seq_len(nrow(sig)) + 1))
rownames(expr) <- rownames(sig)

# These are the production runner's mouse-symbol conversion/preparation calls.
human <- convert_expression_rownames(expr, species = "mouse", target = "human_symbol")
expr_tme <- prepare_tme_expression(human, is_log = FALSE, species = "human", verbose = FALSE)
# references/CIBERSORT/CIBERSORT.R intersects these row names case-sensitively.
overlap <- intersect(rownames(expr_tme), rownames(sig))
cat("Mouse reference rows:", nrow(sig), "\n")
cat("Input genes before conversion:", nrow(expr), "\n")
cat("Prepared human-symbol expression rows:", nrow(expr_tme), "\n")
cat("Native CIBERSORT case-sensitive overlap:", length(overlap), "\n")
cat("Overlapping IDs:", paste(overlap, collapse = ","), "\n")
cat("Case-insensitive overlap (diagnostic only):",
    length(intersect(toupper(rownames(expr_tme)), toupper(rownames(sig)))), "\n")
stopifnot(length(overlap) < nrow(expr_tme) / 10)
cat("REPRODUCED: prepared HGNC input is paired with the MGI reference.\n")
