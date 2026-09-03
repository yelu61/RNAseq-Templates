# Tests for io_utils.R helpers

test_that("detect_count_columns auto-detects count columns", {
  df <- data.frame(
    gene_name = c("A", "B"),
    gene_biotype = c("protein_coding", "lncRNA"),
    Sample1 = c(10, 20),
    Sample2 = c(30, 40),
    stringsAsFactors = FALSE
  )
  cols <- detect_count_columns(df, gene_name_col = "gene_name")
  expect_equal(cols, c("Sample1", "Sample2"))
})

test_that("detect_count_columns respects explicit count_cols", {
  df <- data.frame(
    gene_name = c("A", "B"),
    Sample1 = c(10, 20),
    Sample2 = c(30, 40),
    Sample3 = c(50, 60),
    stringsAsFactors = FALSE
  )
  cols <- detect_count_columns(df, gene_name_col = "gene_name", count_cols = c("Sample1", "Sample3"))
  expect_equal(cols, c("Sample1", "Sample3"))
})

test_that("validate_sample_design catches length mismatch", {
  expect_error(
    validate_sample_design(
      sample_names = c("S1", "S2"),
      groups = c("Ctrl", "Treat"),
      group_levels = c("Ctrl", "Treat"),
      comparisons = list(c("Treat_vs_Ctrl", "Treat", "Ctrl")),
      count_col_names = c("S1", "S2", "S3")
    ),
    "must equal count columns"
  )
})

test_that("validate_sample_design catches missing comparison groups", {
  expect_error(
    validate_sample_design(
      sample_names = c("S1", "S2"),
      groups = c("Ctrl", "Treat"),
      group_levels = c("Ctrl", "Treat"),
      comparisons = list(c("Treat_vs_Ctrl", "Treat", "Missing")),
      count_col_names = c("S1", "S2")
    ),
    "not present in GROUP_LEVELS"
  )
})

test_that("build_count_matrix deduplicates symbols by mean expression", {
  df <- data.frame(
    gene_name = c("A", "A", "B"),
    S1 = c(10, 20, 5),
    S2 = c(30, 40, 6),
    stringsAsFactors = FALSE
  )
  out <- build_count_matrix(df, gene_name_col = "gene_name",
                            count_col_names = c("S1", "S2"),
                            sample_names = c("S1", "S2"))
  expect_equal(rownames(out), c("A", "B"))
  expect_equal(out["A", "S1"], 20)  # row with higher mean kept
})

test_that("build_count_matrix rejects non-integer counts", {
  df <- data.frame(
    gene_name = c("A"),
    S1 = c(10.5),
    stringsAsFactors = FALSE
  )
  expect_error(
    build_count_matrix(df, gene_name_col = "gene_name",
                       count_col_names = "S1", sample_names = "S1"),
    "non-integer"
  )
})

test_that("filter_low_count_genes keeps genes with enough counts", {
  mat <- data.frame(
    S1 = c(10, 0, 5),
    S2 = c(12, 0, 6),
    S3 = c(11, 0, 5),
    S4 = c(0, 0, 5),
    row.names = c("A", "B", "C"),
    stringsAsFactors = FALSE
  )
  out <- filter_low_count_genes(mat, groups = c("G1", "G1", "G2", "G2"), min_count = 5)
  expect_equal(rownames(out$count_data), c("A", "C"))
})

test_that("validate_sample_design stops when shared names are in a different order", {
  expect_error(
    validate_sample_design(
      sample_names = c("S2", "S1"),
      groups = c("Ctrl", "Treat"),
      group_levels = c("Ctrl", "Treat"),
      comparisons = list(c("Treat_vs_Ctrl", "Treat", "Ctrl")),
      count_col_names = c("S1", "S2")
    ),
    "different order"
  )
})

test_that("validate_sample_design accepts matching order of shared names", {
  expect_silent(
    validate_sample_design(
      sample_names = c("S1", "S2"),
      groups = c("Ctrl", "Treat"),
      group_levels = c("Ctrl", "Treat"),
      comparisons = list(c("Treat_vs_Ctrl", "Treat", "Ctrl")),
      count_col_names = c("S1", "S2")
    )
  )
})

test_that("align_to_sample_order reorders named vectors by sample name", {
  aligned <- align_to_sample_order(c(S2 = "b", S1 = "a"), c("S1", "S2"), "GROUPS")
  expect_equal(aligned, c("a", "b"))
  expect_null(names(aligned))
})

test_that("align_to_sample_order passes unnamed vectors through unchanged", {
  expect_equal(align_to_sample_order(c("a", "b"), c("S1", "S2")), c("a", "b"))
})

test_that("align_to_sample_order stops on duplicated or mismatched names", {
  expect_error(
    align_to_sample_order(c(S1 = "a", S3 = "b"), c("S1", "S2"), "BATCH_VECTOR"),
    "do not match the sample names"
  )
  expect_error(
    align_to_sample_order(c(S1 = "a", S1 = "b"), c("S1", "S2")),
    "duplicated names"
  )
})

test_that("build_count_matrix maps named sample_names by column name, not position", {
  df <- data.frame(
    gene_name = c("A", "B"),
    raw_B = c(1, 2),
    raw_A = c(3, 4),
    stringsAsFactors = FALSE
  )
  out <- build_count_matrix(df, gene_name_col = "gene_name",
                            count_col_names = c("raw_B", "raw_A"),
                            sample_names = c(raw_B = "SampleB", raw_A = "SampleA"))
  expect_equal(colnames(out), c("SampleB", "SampleA"))
  expect_equal(out["A", "SampleB"], 1)
  expect_equal(out["A", "SampleA"], 3)
})

test_that("build_count_matrix stops when names(SAMPLE_NAMES) do not match count columns", {
  df <- data.frame(gene_name = "A", S1 = 10, S2 = 20, stringsAsFactors = FALSE)
  expect_error(
    build_count_matrix(df, gene_name_col = "gene_name",
                       count_col_names = c("S1", "S2"),
                       sample_names = c(S1 = "X", WRONG = "Y")),
    "must match the count column names"
  )
})
