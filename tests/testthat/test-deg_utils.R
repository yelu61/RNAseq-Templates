# Tests for deg_utils.R helpers

test_that("validate_threshold_grid accepts valid grid", {
  tg <- data.frame(
    name = c("strict", "loose"),
    p_cutoff = c(0.01, 0.05),
    log2fc = c(1.5, 0.5),
    stringsAsFactors = FALSE
  )
  out <- validate_threshold_grid(tg, "loose")
  expect_equal(out$name, c("strict", "loose"))
})

test_that("validate_threshold_grid rejects duplicated names", {
  tg <- data.frame(
    name = c("strict", "strict"),
    p_cutoff = c(0.01, 0.05),
    log2fc = c(1.5, 0.5),
    stringsAsFactors = FALSE
  )
  expect_error(validate_threshold_grid(tg), "unique")
})

test_that("validate_threshold_grid rejects missing default", {
  tg <- data.frame(
    name = c("strict", "loose"),
    p_cutoff = c(0.01, 0.05),
    log2fc = c(1.5, 0.5),
    stringsAsFactors = FALSE
  )
  expect_error(validate_threshold_grid(tg, "missing"), "DEFAULT_THRESHOLD")
})

test_that("validate_threshold_grid maps legacy padj column", {
  tg <- data.frame(
    name = "standard",
    padj = 0.05,
    log2fc = 1,
    stringsAsFactors = FALSE
  )
  out <- validate_threshold_grid(tg)
  expect_true("p_cutoff" %in% colnames(out))
  expect_equal(out$p_cutoff, 0.05)
})

test_that("mark_deg_by_threshold labels DEGs correctly", {
  res <- data.frame(
    gene_name = c("A", "B", "C", "D", "E"),
    padj = c(0.001, 0.001, 0.5, 0.001, 0.001),
    log2FoldChange = c(2, -2, 2, 0.5, -0.5),
    stringsAsFactors = FALSE
  )
  out <- mark_deg_by_threshold(res, pvalue_thresh = 0.05, log2fc_thresh = 1)
  expect_equal(out$significance, c("Up", "Down", "Not_Sig", "Not_Sig", "Not_Sig"))
})

test_that("genes_for_enrichment splits genes correctly", {
  res <- data.frame(
    gene_name = c("A", "B", "C", "D", "E"),
    padj = c(0.001, 0.001, 0.5, 0.001, 0.001),
    log2FoldChange = c(2, -2, 2, 0.5, -0.5),
    stringsAsFactors = FALSE
  )
  genes <- genes_for_enrichment(res, pvalue_thresh = 0.05, log2fc_thresh = 1)
  expect_equal(sort(genes$sig), c("A", "B"))
  expect_equal(genes$up, "A")
  expect_equal(genes$down, "B")
})

test_that("ranked_gene_list returns sorted named vector", {
  res <- data.frame(
    gene_name = c("A", "B", "C"),
    stat = c(3, 1, 2),
    stringsAsFactors = FALSE
  )
  rnk <- ranked_gene_list(res, rank_column = "stat")
  expect_equal(names(rnk), c("A", "C", "B"))
  expect_equal(unname(rnk), c(3, 2, 1))
})
