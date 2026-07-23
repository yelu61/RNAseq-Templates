# Regenerate the TME deconvolution demo dataset.
# Produces a featureCounts-like table with human gene symbols + gene start/end
# coordinates so counts can be converted to TPM entirely offline (no biomaRt).
# Run from repository root: Rscript examples/demo_RNAseq_TME/regenerate_demo_data.R

set.seed(20250613)

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_file) == 0 || this_file == "") this_file <- "examples/demo_RNAseq_TME/regenerate_demo_data.R"
demo_dir <- dirname(this_file)
data_dir <- file.path(demo_dir, "0-Data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages(library(org.Hs.eg.db))

n_genes <- 6000
samples <- c(paste0("Control_", 1:3), paste0("Treatment_", 1:3))

# Use real human gene symbols so immune signatures and deconvolution reference
# genes have a chance to match.
all_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
genes <- sample(all_symbols, n_genes)

# Ensure key immune-marker genes are present so ssGSEA signatures are non-empty.
marker_genes <- c("CD3D", "CD3E", "CD8A", "CD8B", "CD4", "CD19", "MS4A1", "CD68",
                  "CD14", "NCAM1", "NKG7", "PRF1", "GZMB", "FOXP3", "IL2RA", "CD163",
                  "MRC1", "KIT", "MPO", "S100A8", "S100A9", "C1QA", "CSF1R", "HLA-DRA")
genes <- union(marker_genes, genes)
genes <- genes[seq_len(min(n_genes, length(genes)))]

# Random but valid gene coordinates (start < end), giving kb-scale lengths.
gene_start <- sample(1e6:5e7, length(genes), replace = TRUE)
gene_end   <- gene_start + sample(800:8000, length(genes), replace = TRUE)

# Negative-binomial-like counts; Treatment gets an immune-infiltration shift so
# downstream group comparisons are non-trivial.
base_mu <- 2^runif(length(genes), 3, 12)
counts <- sapply(seq_along(samples), function(j) {
  mu <- base_mu
  if (samples[j] %in% paste0("Treatment_", 1:3)) {
    up <- genes %in% marker_genes
    mu[up] <- mu[up] * 3
  }
  stats::rnbinom(length(genes), mu = mu, size = 4)
})
colnames(counts) <- samples

out <- data.frame(
  gene_name = genes,
  gene_biotype = "protein_coding",
  gene_start = gene_start,
  gene_end = gene_end,
  counts,
  check.names = FALSE
)

write.table(out, file.path(data_dir, "counts.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

meta <- data.frame(
  sample = samples,
  condition = rep(c("Control", "Treatment"), each = 3),
  stringsAsFactors = FALSE
)
write.csv(meta, file.path(data_dir, "metadata.csv"), row.names = FALSE, quote = FALSE)

cat("Wrote", nrow(out), "genes x", length(samples), "samples to", file.path(data_dir, "counts.tsv"), "\n")
cat("Marker genes present:", sum(marker_genes %in% out$gene_name), "/", length(marker_genes), "\n")
