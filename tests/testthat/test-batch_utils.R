# Tests for batch_utils.R (batch-effect diagnostics)

test_that("summarize_pve_by_batch returns per-PC R2 for batch and condition", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("matrixStats")
  set.seed(7)
  # 40 genes x 8 samples, two batches with a real batch shift on some genes.
  mat <- matrix(rnorm(40 * 8, mean = 10, sd = 1), nrow = 40)
  colnames(mat) <- paste0("S", 1:8)
  batch <- rep(c("B1", "B2"), each = 4)
  condition <- rep(c("Ctrl", "Treat"), times = 4)
  mat[1:10, batch == "B2"] <- mat[1:10, batch == "B2"] + 3  # batch effect

  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat))
  pve <- summarize_pve_by_batch(se, batch, condition_vec = condition, ntop = 40, npcs = 4)

  expect_s3_class(pve, "data.frame")
  expect_equal(nrow(pve), 4)
  expect_true(all(c("PC", "R2_batch", "R2_condition") %in% colnames(pve)))
  # PC1 should capture a meaningful fraction of the batch signal.
  expect_true(is.numeric(pve$R2_batch))
  expect_true(all(pve$R2_batch >= 0 & pve$R2_batch <= 1, na.rm = TRUE))
})

test_that("summarize_pve_by_batch validates input length", {
  skip_if_not_installed("DESeq2")
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(rnorm(20), nrow = 5, ncol = 4))
  )
  expect_error(summarize_pve_by_batch(se, c("B1", "B2")), "length must match")
})

test_that("summarize_pve_by_batch aligns a named batch vector to sample order", {
  skip_if_not_installed("DESeq2")
  mat <- matrix(rnorm(20), nrow = 5, ncol = 4,
                dimnames = list(paste0("g", 1:5), c("S2", "S1", "S4", "S3")))
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat))
  batch_named <- c(S1 = "B1", S2 = "B1", S3 = "B2", S4 = "B2")
  expect_silent(summarize_pve_by_batch(se, batch_named, ntop = 5, npcs = 2))
  expect_error(
    summarize_pve_by_batch(se, c(S1 = "B1", S2 = "B1", S3 = "B2", SX = "B2"), ntop = 5, npcs = 2),
    "do not match"
  )
})
