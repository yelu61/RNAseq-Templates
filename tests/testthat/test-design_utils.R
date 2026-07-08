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
