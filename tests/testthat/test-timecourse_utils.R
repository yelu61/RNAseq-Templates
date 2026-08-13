# Tests for timecourse_utils.R (offline helpers; Mfuzz-dependent paths skipped)

test_that("prepare_mfuzz_eset does not open a default graphics device", {
  # Mfuzz imports tcltk, which can emit a harmless XQuartz warning on headless
  # macOS before testthat evaluates the skip condition.
  suppressWarnings(skip_if_not_installed("Mfuzz"))
  skip_if_not_installed("Biobase")
  expr <- matrix(seq_len(24), nrow = 6,
                 dimnames = list(paste0("g", 1:6), paste0("t", 1:4)))
  old_device <- getOption("device")
  on.exit(options(device = old_device), add = TRUE)
  options(device = function(...) stop("unexpected default graphics device"))
  expect_no_error(prepare_mfuzz_eset(expr, min_std = 0))
})

test_that("aggregate_expr_by_group collapses samples by group mean", {
  expr <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE,
                 dimnames = list(c("g1", "g2"), c("s1", "s2", "s3")))
  grp <- c("A", "A", "B")
  out <- aggregate_expr_by_group(expr, grp, fun = mean)
  expect_equal(colnames(out), c("A", "B"))
  expect_equal(out["g1", "A"], mean(c(1, 2)))
  expect_equal(out["g1", "B"], 3)
  expect_equal(out["g2", "A"], mean(c(4, 5)))
  expect_equal(out["g2", "B"], 6)
})

test_that("summarize_timepoint_deg counts up/down DEGs per comparison", {
  res <- data.frame(
    padj = c(0.01, 0.001, 0.2, 0.03, NA),
    log2FoldChange = c(1.0, -1.2, 0.1, 2.0, 5.0),
    stringsAsFactors = FALSE
  )
  res_list <- list(T2_vs_T0 = res, T3_vs_T0 = res)
  out <- summarize_timepoint_deg(res_list, padj_cutoff = 0.05, lfc_cutoff = 0.5)

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2)
  expect_equal(out$Comparison, c("T2_vs_T0", "T3_vs_T0"))
  # rows: (1.0 up), (-1.2 down), (0.1 no), (2.0 up), (NA excluded)
  expect_equal(out$UP[1], 2)
  expect_equal(out$DOWN[1], 1)
  expect_equal(out$Total[1], 3)
})

test_that("summarize_mfuzz_clusters tallies cluster sizes", {
  skip_if_not_installed("dplyr")
  cluster_df <- data.frame(
    gene = paste0("g", 1:6),
    cluster = c(1, 1, 2, 2, 2, 3),
    membership = c(0.9, 0.8, 0.7, 0.6, 0.9, 0.5),
    stringsAsFactors = FALSE
  )
  out <- summarize_mfuzz_clusters(cluster_df)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("cluster", "n_genes") %in% colnames(out)))
  expect_equal(out$n_genes[out$cluster == 1], 2)
  expect_equal(out$n_genes[out$cluster == 2], 3)
  expect_equal(out$n_genes[out$cluster == 3], 1)
})
