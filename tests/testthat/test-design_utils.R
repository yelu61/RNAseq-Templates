# Tests for design_utils.R helpers

test_that("formula_term_labels extracts design terms", {
  expect_equal(formula_term_labels(~ batch + condition), c("batch", "condition"))
})

test_that("validate_batch_design warns when batch is not in design", {
  expect_warning(
    validate_batch_design(c("B1", "B1", "B2", "B2"), ~ condition),
    "does not include"
  )
})

test_that("validate_batch_design accepts batch-adjusted design", {
  expect_silent(validate_batch_design(c("B1", "B1", "B2", "B2"), ~ batch + condition))
})

test_that("validate_batch_design rejects missing batch values", {
  expect_error(validate_batch_design(c("B1", NA, "B2"), ~ batch + condition), "NA or empty")
})

test_that("validate_batch_design checks sample length when provided", {
  expect_error(
    validate_batch_design(c("B1", "B2"), ~ batch + condition, sample_names = c("S1", "S2", "S3")),
    "must equal SAMPLE_NAMES"
  )
})

test_that("add_batch_col aligns a named batch vector to colData rows", {
  colData <- data.frame(condition = factor(c("A", "B")), row.names = c("S1", "S2"))
  out <- add_batch_col(colData, c(S2 = "B2", S1 = "B1"))
  expect_equal(as.character(out$batch), c("B1", "B2"))
})

test_that("add_batch_col stops when a named batch vector does not match samples", {
  colData <- data.frame(condition = factor(c("A", "B")), row.names = c("S1", "S2"))
  expect_error(add_batch_col(colData, c(S1 = "B1", S3 = "B2")), "do not match")
})

test_that("make_paired_col_data aligns named pair_id and groups by sample name", {
  out <- make_paired_col_data(
    sample_names = c("S1", "S2", "S3", "S4"),
    groups = c(S2 = "Treat", S1 = "Ctrl", S4 = "Treat", S3 = "Ctrl"),
    group_levels = c("Ctrl", "Treat"),
    pair_id = c(S3 = "P2", S1 = "P1", S2 = "P1", S4 = "P2")
  )
  expect_equal(as.character(out$condition), c("Ctrl", "Treat", "Ctrl", "Treat"))
  expect_equal(as.character(out$pair_id), c("P1", "P1", "P2", "P2"))
})

test_that("make_paired_col_data stops when a named pair_id does not match samples", {
  expect_error(
    make_paired_col_data(
      sample_names = c("S1", "S2"),
      groups = c("Ctrl", "Treat"),
      group_levels = c("Ctrl", "Treat"),
      pair_id = c(S1 = "P1", SX = "P1")
    ),
    "do not match"
  )
})
