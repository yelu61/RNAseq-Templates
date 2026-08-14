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

test_that("run_go_gsea is reproducible with a fixed seed", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("AnnotationDbi")

  # Fixed ranked ENTREZ list; the function's internal seed governs permutation.
  set.seed(42)
  ids <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "ENTREZID")
  ids <- sample(ids, 2000)
  ranked <- sort(stats::setNames(stats::rnorm(length(ids)), ids), decreasing = TRUE)

  r1 <- run_go_gsea(ranked, org_db = org.Hs.eg.db::org.Hs.eg.db, ont = "BP", p_cutoff = 0.05)
  r2 <- run_go_gsea(ranked, org_db = org.Hs.eg.db::org.Hs.eg.db, ont = "BP", p_cutoff = 0.05)
  skip_if(is.null(r1) || is.null(r2), "gseGO unavailable in this environment")
  expect_identical(as.data.frame(r1), as.data.frame(r2))
})

test_that("run_go_gsea returns NULL for a too-short ranked list", {
  skip_if_not_installed("clusterProfiler")
  out <- run_go_gsea(c("1" = 1, "2" = 2), org_db = NULL, min_size = 10)
  expect_null(out)
})
