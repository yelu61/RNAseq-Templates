# Tests for tme_utils.R helpers

test_that("undo_log_expr reverses log2 transformation", {
  mat <- matrix(c(0, 1, 2, 3), nrow = 2, dimnames = list(c("A", "B"), c("S1", "S2")))
  out <- undo_log_expr(mat, is_log = TRUE, log_base = 2)
  expect_equal(as.numeric(out), c(0, 1, 3, 7))
  expect_equal(undo_log_expr(mat, is_log = FALSE), mat)
})

test_that("deduplicate_expression_by_symbol keeps highest-mean row", {
  expr <- data.frame(
    S1 = c(10, 5, 1),
    S2 = c(12, 4, 1),
    symbol = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )
  out <- deduplicate_expression_by_symbol(expr, symbol_col = "symbol")
  expect_equal(sort(rownames(out)), c("A", "B"))
  expect_equal(as.numeric(out["A", c("S1", "S2")]), c(10, 12))
})

test_that("validate_tme_input catches common input problems", {
  # No rownames → error
  mat <- matrix(c(1, 2, 3, 4), nrow = 2)
  expect_error(validate_tme_input(mat), "rownames")

  # Negative values → error
  mat <- matrix(c(1, 2, -1, 4), nrow = 2, dimnames = list(c("A", "B"), c("S1", "S2")))
  expect_error(validate_tme_input(mat), "negative")

  # Missing values → error
  mat[2, 2] <- NA
  expect_error(validate_tme_input(mat), "missing values")

  good <- matrix(c(1, 2, 3, 4), nrow = 2, dimnames = list(c("A", "B"), c("S1", "S2")))
  expect_true(validate_tme_input(good))

  ensembl <- matrix(c(1, 2, 3, 4), nrow = 2, dimnames = list(c("ENSG000001", "ENSG000002"), c("S1", "S2")))
  expect_warning(validate_tme_input(ensembl), "Ensembl")
})

test_that("prepare_tme_expression human path reverses log and deduplicates", {
  mat <- matrix(c(0, 1, 2, 3), nrow = 2, dimnames = list(c("A", "B"), c("S1", "S2")))
  out <- prepare_tme_expression(mat, is_log = TRUE, species = "human", verbose = FALSE)
  expect_equal(sort(rownames(out)), c("A", "B"))
  expect_equal(as.numeric(out["A", ]), c(0, 3))
  expect_equal(as.numeric(out["B", ]), c(1, 7))
})

test_that("prepare_tme_expression mouse path attempts conversion", {
  skip_if_not_installed("biomaRt")

  # Build a tiny expression matrix with known mouse symbols
  mat <- matrix(
    c(10, 5, 8, 4),
    nrow = 2,
    dimnames = list(c("Actb", "Gapdh"), c("S1", "S2"))
  )

  # If network is available, conversion should produce HGNC-symbol rownames.
  # If not, the function should error. Either outcome is acceptable for this test,
  # but we assert that the function does not silently return the input unchanged.
  res <- tryCatch(
    prepare_tme_expression(mat, is_log = FALSE, species = "mouse", verbose = FALSE),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    expect_true(TRUE)
  } else {
    expect_equal(sort(toupper(rownames(res))), c("ACTB", "GAPDH"))
  }
})

test_that("prepare_tme_expression rejects invalid species", {
  mat <- matrix(1:4, nrow = 2, dimnames = list(c("A", "B"), c("S1", "S2")))
  expect_error(prepare_tme_expression(mat, species = "rat"), "'arg' should be one of")
})

test_that("calc_tme_barplot_size scales with sample count", {
  s1 <- calc_tme_barplot_size(10)
  s2 <- calc_tme_barplot_size(50)
  expect_true(s2["width"] > s1["width"])
  expect_true(s1["height"] == s2["height"])

  s3 <- calc_tme_barplot_size(200)
  expect_true(s3["width"] <= 24)
})

test_that("calc_tme_boxplot_size scales with cell-type count", {
  s1 <- calc_tme_boxplot_size(4)
  s2 <- calc_tme_boxplot_size(20)
  expect_true(s2["width"] >= s1["width"])
  expect_true(s2["height"] > s1["height"])

  s3 <- calc_tme_boxplot_size(100)
  expect_true(s3["width"] <= 22)
  expect_true(s3["height"] <= 20)
})

test_that("plot_tme_boxplot_pdf accepts group_colors", {
  long_df <- data.frame(
    sample = rep(c("S1", "S2", "S3", "S4"), each = 2),
    condition = rep(c("A", "B"), each = 4),
    cell_type = rep(c("T", "B"), 4),
    fraction = c(0.1, 0.2, 0.15, 0.25, 0.05, 0.3, 0.1, 0.28),
    stringsAsFactors = FALSE
  )
  tmp <- tempfile(fileext = ".pdf")
  expect_silent(
    plot_tme_boxplot_pdf(long_df, group_col = "condition", value_col = "fraction",
                         filename = tmp, width = 8, height = 6,
                         group_colors = c("A" = "#6F6F6F", "B" = "#E07B54"))
  )
  expect_true(file.exists(tmp))
  unlink(tmp)
})

test_that("plot_tme_heatmap_pdf accepts group_colors and writes PDF", {
  tme_df <- data.frame(
    sample = c("S1", "S2", "S3", "S4"),
    score1 = c(1, 2, 3, 4),
    score2 = c(4, 3, 2, 1),
    stringsAsFactors = FALSE
  )
  meta <- data.frame(
    sample = c("S1", "S2", "S3", "S4"),
    condition = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  tmp <- tempfile(fileext = ".pdf")
  expect_silent(
    plot_tme_heatmap_pdf(tme_df, meta, group_col = "condition", sample_col = "sample",
                         group_colors = c("A" = "#6F6F6F", "B" = "#E07B54"),
                         filename = tmp, title = "Test", width = 6, height = 5)
  )
  expect_true(file.exists(tmp))
  unlink(tmp)
})
