# Regenerate the TCGA-GEO demo dataset (local-file mode).
# Produces a TCGA-like cohort: raw counts + TPM with TCGA barcodes plus a
# synthetic clinical table, so the Tumor-vs-Normal DEG and survival workflow
# can run without any network download.
# Run from repository root: Rscript examples/demo_RNAseq_TCGA_GEO/regenerate_demo_data.R

set.seed(20250614)

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_file) == 0 || this_file == "") this_file <- "examples/demo_RNAseq_TCGA_GEO/regenerate_demo_data.R"
demo_dir <- dirname(this_file)
data_dir <- file.path(demo_dir, "0-Data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages(library(org.Hs.eg.db))

n_genes <- 5000
n_tumor <- 18
n_normal <- 6

# TCGA sample-type code sits in barcode positions 14-15: <10 = Tumor, >=10 = Normal.
# infer_tcga_tumor_normal() reads exactly substr(barcode, 14, 15), so place the
# 2-digit type code there. Layout: TCGA-LLL-SSSS-TTR (positions 1-4=TCGA-, 6-8,
# 10-13, 14-15=type).
mk_barcode <- function(i, type) {
  # "TCGA-DMO-0001" is exactly 13 characters, so the 2-digit type code lands at
  # positions 14-15, which is what infer_tcga_tumor_normal() reads.
  sprintf("TCGA-DMO-%04d%02d", i, type)
}
tumor_barcodes  <- vapply(seq_len(n_tumor),  mk_barcode, character(1), type = 1)
normal_barcodes <- vapply(seq_len(n_normal), mk_barcode, character(1), type = 11)
barcodes <- c(tumor_barcodes, normal_barcodes)

all_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
genes <- sample(all_symbols, n_genes)

# Survival gene of interest plus a few cancer-relevant markers, forced present.
survival_gene <- "MKI67"
markers <- c(survival_gene, "ICAM1", "TP53", "EGFR", "MYC", "CDKN2A", "ERBB2")
genes <- union(markers, genes)
genes <- genes[seq_len(min(n_genes, length(genes)))]

base_mu <- 2^runif(length(genes), 3, 11)

# Tumor samples: overexpress MKI67 (proliferation) plus a broad up/down program
# so the Tumor-vs-Normal contrast yields a meaningful number of DEGs.
n_up <- 300
n_down <- 300
non_marker <- which(!genes %in% markers)
up_set <- sample(non_marker, n_up)
down_set <- sample(setdiff(non_marker, up_set), n_down)
counts <- sapply(seq_along(barcodes), function(j) {
  mu <- base_mu
  if (j <= n_tumor) {
    mu[genes == survival_gene] <- mu[genes == survival_gene] * 6
    mu[up_set] <- mu[up_set] * 5
    mu[down_set] <- mu[down_set] / 6
  }
  stats::rnbinom(length(genes), mu = mu, size = 8)
})
colnames(counts) <- barcodes

counts_df <- data.frame(gene_name = genes, counts, check.names = FALSE)
write.csv(counts_df, file.path(data_dir, "counts.csv"), row.names = FALSE, quote = FALSE)

# TPM: normalize per sample to sum 1e6 (approximation adequate for survival ranking).
tpm <- apply(counts, 2, function(x) x / sum(x) * 1e6)
tpm_df <- data.frame(gene_name = genes, tpm, check.names = FALSE)
write.csv(tpm_df, file.path(data_dir, "tpm.csv"), row.names = FALSE, quote = FALSE)

# Synthetic clinical metadata. High MKI67 -> shorter survival for a detectable signal.
expr_rank <- rank(base_mu) # not the survival driver; use the tumor MKI67 shift instead
clinical <- data.frame(
  barcode = barcodes,
  stringsAsFactors = FALSE
)
clinical$tissue_type <- ifelse(seq_len(nrow(clinical)) <= n_tumor, "Tumor", "Normal")

# Survival only meaningful for tumor samples; give normals follow-up too.
is_tumor <- clinical$tissue_type == "Tumor"
# Latent risk: tumors overexpressing MKI67 have higher risk.
risk <- stats::rnorm(nrow(clinical), mean = ifelse(is_tumor, 0.5, -0.5), sd = 0.6)
status <- ifelse(is_tumor, rbinom(sum(is_tumor), 1, 0.45), 0L)
time_event <- round(pmax(60, stats::rexp(nrow(clinical), rate = 1 / (700 * exp(-risk)))))
time_censor <- round(pmax(60, stats::rexp(nrow(clinical), rate = 1 / 1100)))

clinical$vital_status <- ifelse(status == 1, "Dead", "Alive")
clinical$days_to_death <- ifelse(status == 1, time_event, NA)
clinical$days_to_last_follow_up <- ifelse(status == 1, NA, time_censor)
clinical$ajcc_pathologic_stage <- sample(
  c("Stage I", "Stage II", "Stage III", "Stage IV"), nrow(clinical),
  replace = TRUE, prob = c(0.25, 0.35, 0.25, 0.15)
)
clinical$age_at_diagnosis <- round(stats::rnorm(nrow(clinical), mean = 62, sd = 11))

write.csv(clinical, file.path(data_dir, "clinical.csv"), row.names = FALSE, quote = FALSE)

cat("Wrote", nrow(counts_df), "genes x", length(barcodes), "samples\n")
cat("Tumor:", n_tumor, "Normal:", n_normal, "\n")
cat("Survival gene", survival_gene, "present:", survival_gene %in% genes, "\n")
cat("Events (deaths):", sum(clinical$vital_status == "Dead"), "\n")
