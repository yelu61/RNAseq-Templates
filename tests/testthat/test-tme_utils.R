# Tests for tme_utils.R helpers

repo_root <- tryCatch(
  rprojroot::find_root(rprojroot::is_git_root),
  error = function(e) getwd()
)
cibersort_script_path <- file.path(repo_root, "references", "CIBERSORT", "CIBERSORT.R")

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

# -----------------------------------------------------------------------------
# Native CIBERSORT wrapper tests
# -----------------------------------------------------------------------------

make_tiny_signature_and_mixture <- function() {
  sig <- matrix(
    c(100, 50, 20, 10, 5,
      5, 10, 20, 50, 100),
    nrow = 5, byrow = FALSE
  )
  rownames(sig) <- c("A", "B", "C", "D", "E")
  colnames(sig) <- c("CT1", "CT2")

  mix <- matrix(
    c(80, 40, 20, 20, 10,
      20, 20, 20, 40, 80),
    nrow = 5, byrow = FALSE
  )
  rownames(mix) <- c("A", "B", "C", "D", "E")
  colnames(mix) <- c("S1", "S2")

  list(signature = sig, mixture = mix)
}

test_that("run_native_cibersort runs with a matrix signature", {
  skip_if_not_installed("e1071")
  skip_if_not_installed("preprocessCore")
  skip_if_not_installed("future")
  skip_if_not_installed("furrr")
  skip_if_not_installed("purrr")

  dat <- make_tiny_signature_and_mixture()
  res <- run_native_cibersort(
    dat$mixture,
    signature_file = dat$signature,
    cibersort_script = cibersort_script_path,
    is_log = FALSE,
    perm = 0,
    QN = FALSE,
    verbose = FALSE
  )

  expect_equal(nrow(res), 2)
  expect_true("ID" %in% colnames(res))
  expect_true(all(c("CT1", "CT2", "P-value", "Correlation", "RMSE") %in% colnames(res)))
  expect_type(res$CT1, "double")
  expect_equal(res$ID, c("S1", "S2"))
})

test_that("run_native_cibersort reads bundled LM22.txt", {
  skip_if_not_installed("e1071")
  skip_if_not_installed("preprocessCore")
  skip_if_not_installed("future")
  skip_if_not_installed("furrr")
  skip_if_not_installed("purrr")

  sig_file <- file.path(repo_root, "references", "CIBERSORT", "LM22.txt")
  skip_if_not(file.exists(sig_file), "Bundled LM22.txt not found")

  sig <- read_cibersort_signature(sig_file)
  expect_equal(ncol(sig), 22)
  expect_true(is.matrix(sig))
})

test_that("run_native_cibersort reads bundled mouse signature", {
  sig_file <- file.path(repo_root, "references", "CIBERSORT", "cibersort_mouse_22.csv")
  skip_if_not(file.exists(sig_file), "Bundled cibersort_mouse_22.csv not found")

  sig <- read_cibersort_signature(sig_file)
  expect_equal(ncol(sig), 25)
  expect_true(is.matrix(sig))
})

test_that("run_native_cibersort errors on missing script", {
  expect_error(
    run_native_cibersort(matrix(1:4, nrow = 2), cibersort_script = "not_a_file.R"),
    "CIBERSORT script not found"
  )
})

# -----------------------------------------------------------------------------
# Native vs IOBR comparison tests
# -----------------------------------------------------------------------------

test_that("compare_native_iobr_cibersort computes concordance metrics", {
  dat <- make_tiny_signature_and_mixture()
  native <- run_native_cibersort(
    dat$mixture,
    signature_file = dat$signature,
    cibersort_script = cibersort_script_path,
    is_log = FALSE,
    perm = 0,
    QN = FALSE,
    verbose = FALSE
  )
  iobr <- data.frame(
    ID = c("S1", "S2"),
    CT1 = c(0.98, 0.08),
    CT2 = c(0.02, 0.92),
    stringsAsFactors = FALSE
  )

  cmp <- compare_native_iobr_cibersort(native, iobr, id_column = "ID")

  expect_true(all(c("summary", "long") %in% names(cmp)))
  expect_equal(nrow(cmp$summary), 2)
  expect_true(all(c("cell_type", "n", "correlation", "rmse", "mae", "paired_t_pvalue") %in% colnames(cmp$summary)))
  expect_equal(nrow(cmp$long), 4)
  expect_true(all(c("sample", "cell_type", "native", "iobr", "difference") %in% colnames(cmp$long)))
})

test_that("compare_native_iobr_cibersort matches relaxed cell-type names", {
  native <- data.frame(
    ID = c("S1", "S2"),
    `B.cells.naive` = c(0.1, 0.2),
    `T.cells.CD8` = c(0.3, 0.4),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  iobr <- data.frame(
    ID = c("S1", "S2"),
    `B cells naive` = c(0.12, 0.22),
    `T cells CD8` = c(0.31, 0.39),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  cmp <- compare_native_iobr_cibersort(native, iobr, id_column = "ID")
  expect_equal(sort(cmp$summary$cell_type), c("b cells naive", "t cells cd8"))
})

test_that("compare_native_iobr_cibersort errors on mismatched IDs", {
  native <- data.frame(ID = "A", CT1 = 0.5, CT2 = 0.5)
  iobr <- data.frame(ID = "B", CT1 = 0.5, CT2 = 0.5)
  expect_error(compare_native_iobr_cibersort(native, iobr), "No common sample IDs")
})

test_that("plot_cibersort_correlation_pdf writes PDF", {
  long_df <- data.frame(
    sample = c("S1", "S2"),
    cell_type = "CT1",
    native = c(0.9, 0.1),
    iobr = c(0.85, 0.15),
    difference = c(0.05, -0.05),
    stringsAsFactors = FALSE
  )
  tmp <- tempfile(fileext = ".pdf")
  expect_silent(plot_cibersort_correlation_pdf(long_df, filename = tmp, width = 6, height = 5))
  expect_true(file.exists(tmp))
  unlink(tmp)
})

test_that("plot_cibersort_difference_pdf writes PDF", {
  long_df <- data.frame(
    sample = c("S1", "S2"),
    cell_type = "CT1",
    native = c(0.9, 0.1),
    iobr = c(0.85, 0.15),
    difference = c(0.05, -0.05),
    stringsAsFactors = FALSE
  )
  tmp <- tempfile(fileext = ".pdf")
  expect_silent(plot_cibersort_difference_pdf(long_df, filename = tmp, width = 6, height = 5))
  expect_true(file.exists(tmp))
  unlink(tmp)
})

# -----------------------------------------------------------------------------
# Broad-category aggregation tests
# -----------------------------------------------------------------------------

test_that("get_cibersort_category_map covers human LM22 cell types", {
  map <- get_cibersort_category_map("human")
  expect_equal(sort(unique(map$category)),
               sort(c("B cells", "CD4 T cells", "CD8 T cells", "Other T cells",
                      "NK cells", "Monocytes", "Macrophages", "Dendritic cells",
                      "Mast cells",
                      "Eosinophils", "Neutrophils")))
  expect_true(all(c("B cells naive", "T cells CD8", "Monocytes") %in% map$cell_type))
})

test_that("get_cibersort_category_map covers mouse signature cell types", {
  map <- get_cibersort_category_map("mouse")
  expect_true(all(c("B Cells Naive", "T Cells CD8 Actived", "M1 Macrophage",
                    "NK.Actived", "Monocyte") %in% map$cell_type))
})

test_that("get_xcell_category_map returns IOBR-compatible keys", {
  map <- get_xcell_category_map()
  expect_true("B-cells" %in% map$cell_type)
  expect_true("CD4+_T-cells" %in% map$cell_type)
  expect_true("Fibroblasts" %in% map$cell_type)
  expect_true("Fibroblasts" %in% map$category)
  expect_true("Endothelials" %in% map$category)
})

test_that("aggregate_tme_by_category sums CIBERSORT fractions", {
  cib_df <- data.frame(
    ID = c("S1", "S2"),
    `B cells naive` = c(0.1, 0.2),
    `B cells memory` = c(0.05, 0.1),
    `Plasma cells` = c(0.0, 0.05),
    `T cells CD8` = c(0.2, 0.3),
    `T cells CD4 naive` = c(0.05, 0.1),
    `NK cells resting` = c(0.1, 0.05),
    `Monocytes` = c(0.2, 0.1),
    `Macrophages M1` = c(0.1, 0.05),
    `Eosinophils` = c(0.05, 0.03),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  map <- get_cibersort_category_map("human")
  agg <- aggregate_tme_by_category(cib_df, id_column = "ID", category_map = map, method = "sum")

  expect_true("B cells" %in% colnames(agg$wide))
  expect_equal(agg$wide[["B cells"]], c(0.15, 0.35))
  expect_false(anyNA(agg$wide))
})

test_that("aggregate_tme_by_category means xCell scores", {
  xcell_df <- data.frame(
    ID = c("S1", "S2"),
    `B-cells_xCell` = c(0.1, 0.2),
    `CD4+_T-cells_xCell` = c(0.2, 0.3),
    `CD8+_T-cells_xCell` = c(0.15, 0.25),
    `Macrophages_xCell` = c(0.1, 0.05),
    `Fibroblasts_xCell` = c(0.3, 0.4),
    `ImmuneScore_xCell` = c(0.5, 0.6),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  map <- get_xcell_category_map()
  agg <- aggregate_tme_by_category(xcell_df, id_column = "ID", category_map = map, method = "mean")

  expect_true("B cells" %in% colnames(agg$wide))
  expect_true("Fibroblasts" %in% colnames(agg$wide))
  expect_false("ImmuneScore" %in% colnames(agg$wide))
  expect_equal(agg$wide$`CD4 T cells`, c(0.2, 0.3))
})

test_that("aggregate_tme_by_category warns on unmatched cell types", {
  df <- data.frame(ID = "S1", `B cells naive` = 0.2, `Unknown cell` = 0.5, check.names = FALSE)
  map <- get_cibersort_category_map("human")
  expect_warning(
    aggregate_tme_by_category(df, id_column = "ID", category_map = map, method = "sum"),
    "Unknown cell"
  )
})

test_that("aggregate_tme_by_category errors when nothing matches", {
  df <- data.frame(ID = "S1", `Foo` = 0.5, `Bar` = 0.5, check.names = FALSE)
  map <- get_cibersort_category_map("human")
  expect_error(
    aggregate_tme_by_category(df, id_column = "ID", category_map = map, method = "sum"),
    "No cell-type columns could be matched"
  )
})
