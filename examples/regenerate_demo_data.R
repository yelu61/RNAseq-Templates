# Regenerate demo_counts.tsv with real mouse gene symbols.
# Run from repository root: Rscript examples/regenerate_demo_data.R

set.seed(42)

input_file  <- "examples/demo_data/demo_counts.tsv"
output_file <- "examples/demo_data/demo_counts.tsv"

stopifnot(file.exists(input_file))

raw <- read.delim(input_file, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot("gene_name" %in% colnames(raw))

library(org.Mm.eg.db)
all_symbols <- keys(org.Mm.eg.db, keytype = "SYMBOL")
new_symbols <- sample(all_symbols, nrow(raw))

# Replace gene_name column but keep everything else intact
raw$gene_name <- new_symbols

# Write back as tab-separated with original header order preserved
write.table(
  raw,
  file = output_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Regenerated", output_file, "with", length(unique(new_symbols)), "unique real mouse symbols\n")
