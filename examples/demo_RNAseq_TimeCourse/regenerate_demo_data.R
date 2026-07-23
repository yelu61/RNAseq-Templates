# Regenerate the TimeCourse demo dataset.
# Produces a VST-scale expression matrix plus metadata across four time points,
# and a raw-count table + count metadata for the time-point-vs-baseline DEG.
# Temporal patterns are injected so Mfuzz finds real clusters.
# Run from repository root: Rscript examples/demo_RNAseq_TimeCourse/regenerate_demo_data.R

set.seed(20250615)

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_file) == 0 || this_file == "") this_file <- "examples/demo_RNAseq_TimeCourse/regenerate_demo_data.R"
demo_dir <- dirname(this_file)
data_dir <- file.path(demo_dir, "0-Data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages(library(org.Hs.eg.db))

time_points <- c("Day0", "Day7", "Day14", "Day21")
n_rep <- 3
n_genes <- 3000
n_pattern <- 600  # genes per temporal pattern; the rest are flat

all_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
genes <- sample(all_symbols, n_genes)

# Sample names encode time point and replicate.
samples <- as.vector(outer(time_points, seq_len(n_rep), function(t, r) paste0(t, "_Rep", r)))
meta <- data.frame(
  sample = samples,
  time = factor(rep(time_points, each = n_rep), levels = time_points),
  condition = "Treatment",
  stringsAsFactors = FALSE
)

# Baseline log2-scale expression per gene.
base <- runif(n_genes, 5, 11)

# Temporal patterns (log2-fold change relative to Day0 across the 4 time points).
patterns <- list(
  up        = c(0, 1.0, 2.0, 3.0),
  down      = c(0, -1.0, -2.0, -3.0),
  transient = c(0, 2.5, 1.0, 0),
  late_up   = c(0, 0, 0.5, 3.0)
)
pattern_of <- rep("flat", n_genes)
idx <- seq_len(min(n_pattern * length(patterns), n_genes))
pattern_of[idx] <- rep(names(patterns), each = n_pattern)[seq_along(idx)]

# Build VST-scale expression (no negatives so the scale detector reads it as
# normalized expression, not raw counts or integer counts).
expr <- sapply(seq_along(samples), function(j) {
  tp <- meta$time[j]
  tp_idx <- match(tp, time_points)
  vals <- base + stats::rnorm(n_genes, 0, 0.25)
  for (pn in names(patterns)) {
    sel <- pattern_of == pn
    vals[sel] <- vals[sel] + patterns[[pn]][tp_idx]
  }
  pmax(vals, 0.5)
})
colnames(expr) <- samples
rownames(expr) <- genes

expr_df <- data.frame(gene_name = genes, expr, check.names = FALSE)
write.csv(expr_df, file.path(data_dir, "vsd_matrix.csv"), row.names = FALSE, quote = FALSE)
write.csv(meta, file.path(data_dir, "colData.csv"), row.names = FALSE, quote = FALSE)

# Raw integer counts consistent with the same temporal structure, for DEG.
counts <- sapply(seq_along(samples), function(j) {
  tp <- meta$time[j]
  tp_idx <- match(tp, time_points)
  log_mu <- base + stats::rnorm(n_genes, 0, 0.15)
  for (pn in names(patterns)) {
    sel <- pattern_of == pn
    log_mu[sel] <- log_mu[sel] + patterns[[pn]][tp_idx]
  }
  mu <- 2^log_mu
  stats::rnbinom(n_genes, mu = mu, size = 6)
})
colnames(counts) <- samples
counts_df <- data.frame(gene_name = genes, gene_biotype = "protein_coding", counts, check.names = FALSE)
write.table(counts_df, file.path(data_dir, "raw_counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(meta, file.path(data_dir, "metadata.csv"), row.names = FALSE, quote = FALSE)

cat("Wrote", nrow(expr_df), "genes x", length(samples), "samples across", length(time_points), "time points\n")
cat("Pattern genes:", sum(pattern_of != "flat"), "| flat:", sum(pattern_of == "flat"), "\n")
