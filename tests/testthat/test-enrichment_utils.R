# Tests for enrichment_utils.R helpers

test_that("enrich_result_to_df returns empty df for NULL", {
  out <- enrich_result_to_df(NULL)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0)
})

test_that("enrich_result_to_df returns empty df for empty object", {
  out <- enrich_result_to_df(data.frame())
  expect_equal(nrow(out), 0)
})

test_that("build_multi_comparison_enrich_df handles empty list", {
  out <- build_multi_comparison_enrich_df(list())
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0)
})

test_that("build_multi_comparison_enrich_df binds comparison column", {
  df1 <- data.frame(
    ID = "GO:001",
    Description = "term1",
    p.adjust = 0.01,
    Count = 5,
    stringsAsFactors = FALSE
  )
  df2 <- data.frame(
    ID = "GO:002",
    Description = "term2",
    p.adjust = 0.02,
    Count = 3,
    stringsAsFactors = FALSE
  )
  out <- build_multi_comparison_enrich_df(list(A = df1, B = df2))
  expect_equal(nrow(out), 2)
  expect_equal(out$comparison, c("A", "B"))
})

test_that("build_multi_comparison_enrich_df applies ontology filter", {
  df1 <- data.frame(
    ID = "GO:001",
    Description = "term1",
    p.adjust = 0.01,
    Count = 5,
    ONTOLOGY = "BP",
    stringsAsFactors = FALSE
  )
  df2 <- data.frame(
    ID = "GO:002",
    Description = "term2",
    p.adjust = 0.02,
    Count = 3,
    ONTOLOGY = "MF",
    stringsAsFactors = FALSE
  )
  out <- build_multi_comparison_enrich_df(list(A = df1, B = df2), ontology_filter = "BP")
  expect_equal(nrow(out), 1)
  expect_equal(out$ID, "GO:001")
})
