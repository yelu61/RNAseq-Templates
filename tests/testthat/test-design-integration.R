# Integration tests for the General paired-design and batch-correction paths.
#
# These reproduce the exact failure modes fixed in run_analysis.R / RNAseq_General.ipynb:
#   1. PAIR_ID set -> downstream filter_low_count_genes()/save() need `group`, which the
#      paired branch previously never created ("object 'group' not found").
#   2. BATCH_VECTOR set with DESIGN_FORMULA = ~ batch + condition -> colData must actually
#      contain a `batch` column, otherwise DESeqDataSetFromMatrix fails to resolve the term.

# Build a small synthetic count matrix with sample-named columns.
make_test_counts <- function(n_genes = 200, sample_names) {
  set.seed(42)
  m <- matrix(
    rpois(n_genes * length(sample_names), lambda = 50),
    nrow = n_genes,
    dimnames = list(paste0("gene", seq_len(n_genes)), sample_names)
  )
  as.data.frame(m)
}

test_that("add_batch_col appends a batch factor and validates length", {
  samples <- paste0("S", 1:4)
  cd <- data.frame(row.names = samples, sample = samples,
                   condition = factor(rep(c("A", "B"), each = 2)))
  out <- add_batch_col(cd, c("B1", "B1", "B2", "B2"))
  expect_true("batch" %in% colnames(out))
  expect_s3_class(out$batch, "factor")
  expect_equal(as.character(out$batch), c("B1", "B1", "B2", "B2"))
  expect_error(add_batch_col(cd, c("B1", "B2")), "does not match")
  expect_identical(add_batch_col(cd, NULL), cd)
})

test_that("paired branch yields group for downstream filtering", {
  samples <- paste0(rep(c("Ctrl", "Treat"), each = 3), "_", 1:3)
  groups <- rep(c("Ctrl", "Treat"), each = 3)
  pair_id <- rep(1:3, times = 2)
  group_levels <- c("Ctrl", "Treat")

  countData <- make_test_counts(sample_names = samples)

  # Mirror the fixed paired branch logic.
  validate_paired_design(samples, groups, pair_id, group_levels)
  colData <- make_paired_col_data(samples, groups, group_levels, pair_id)
  group <- colData$condition

  # This is the call that previously failed with object 'group' not found.
  expect_error(filter_res <- filter_low_count_genes(countData, group, 10), NA)
  expect_true(nrow(filter_res$count_data) > 0)
  expect_equal(nrow(colData), ncol(countData))
})

test_that("paired + batch produces a colData that resolves the design formula", {
  samples <- paste0(rep(c("Ctrl", "Treat"), each = 4), "_", 1:4)
  groups <- rep(c("Ctrl", "Treat"), each = 4)
  pair_id <- rep(1:4, times = 2)
  batch <- rep(c("B1", "B2"), length.out = length(samples))
  group_levels <- c("Ctrl", "Treat")

  countData <- make_test_counts(sample_names = samples)

  validate_paired_design(samples, groups, pair_id, group_levels)
  colData <- make_paired_col_data(samples, groups, group_levels, pair_id)
  colData <- add_batch_col(colData, batch)

  # batch + condition formula must resolve against colData columns.
  design <- ~ batch + condition
  expect_true(all(formula_term_labels(design) %in% colnames(colData)))
  # model matrix construction should succeed for every sample.
  mm <- stats::model.matrix(design, data = colData)
  expect_equal(nrow(mm), length(samples))
})

test_that("non-paired batch path writes batch into colData", {
  samples <- paste0("S", 1:6)
  groups <- rep(c("A", "B"), each = 3)
  group_levels <- c("A", "B")
  batch <- rep(c("B1", "B2", "B3"), times = 2)

  countData <- make_test_counts(sample_names = samples)
  colData <- make_col_data(countData, samples, groups, group_levels)
  colData <- add_batch_col(colData, batch)

  expect_true(all(c("condition", "batch") %in% colnames(colData)))
  mm <- stats::model.matrix(~ batch + condition, data = colData)
  expect_equal(nrow(mm), length(samples))
})
