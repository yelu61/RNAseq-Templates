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

test_that("plot_deg_upset_pdf renders a non-empty PDF in a non-interactive closure", {
  testthat::skip_if_not_installed("UpSetR")
  # Regression: UpSetR::upset() returns a grid object that is only auto-printed at
  # the top level; inside the save closure it must be print()ed or the device stays
  # blank and .promote_validated_figure() refuses the <1500-byte PDF.
  sets <- list(
    S1 = c("g1", "g2", "g3", "g9"),
    S2 = c("g2", "g3", "g4"),
    S3 = c("g3", "g5", "g6")
  )
  outfile <- file.path(tempdir(), "test_deg_upset.pdf")
  on.exit(unlink(outfile), add = TRUE)
  expect_invisible(plot_deg_upset_pdf(sets, outfile, title = "test"))
  expect_true(file.exists(outfile))
  expect_gt(file.info(outfile)$size, 1500)
})

test_that("validate_threshold_grid validates per-row p_column", {
  tg <- data.frame(
    name = c("standard", "exploratory"),
    p_cutoff = c(0.05, 0.05),
    log2fc = c(1.0, 0.5),
    p_column = c("padj", "pvalue"),
    stringsAsFactors = FALSE
  )
  expect_silent(validate_threshold_grid(tg, "standard"))
  tg$p_column[2] <- "fdr"
  expect_error(validate_threshold_grid(tg), "p_column")
})

test_that("threshold_p_column resolves per-row override with default fallback", {
  th <- data.frame(name = "exploratory", p_cutoff = 0.05, log2fc = 0.5,
                   p_column = "pvalue", stringsAsFactors = FALSE)
  expect_equal(threshold_p_column(th, "padj"), "pvalue")
  th$p_column <- NA
  expect_equal(threshold_p_column(th, "padj"), "padj")
  th$p_column <- NULL
  expect_equal(threshold_p_column(th, "padj"), "padj")
  expect_equal(threshold_p_column(th, NULL), "padj")
})

test_that("build_deg_threshold_sets honors a mixed padj/pvalue grid", {
  res <- data.frame(
    gene_name = c("A", "B", "C", "D"),
    pvalue = c(0.001, 0.01, 0.03, 0.20),
    padj   = c(0.04, 0.30, 0.60, 0.90),
    log2FoldChange = c(2, 1.5, -1.2, 0.8),
    stringsAsFactors = FALSE
  )
  grid <- data.frame(
    name = c("standard", "exploratory"),
    p_cutoff = c(0.05, 0.05),
    log2fc = c(1, 1),
    p_column = c("padj", "pvalue"),
    stringsAsFactors = FALSE
  )
  sets <- build_deg_threshold_sets(list(T_vs_C = res), grid, pvalue_column = "padj")
  # padj tier: only gene A passes; nominal-p tier: A, B, C pass.
  expect_equal(sum(sets$standard$T_vs_C$significance != "Not_Sig"), 1)
  expect_equal(sum(sets$exploratory$T_vs_C$significance != "Not_Sig"), 3)
  # Summary records the per-row column and definition.
  s <- summarize_deg_thresholds(sets, grid, pvalue_column = "padj", lfc_column = "log2FoldChange")
  expect_equal(s$pvalue_column[s$Threshold == "standard"], "padj")
  expect_equal(s$pvalue_column[s$Threshold == "exploratory"], "pvalue")
  expect_match(s$definition[s$Threshold == "exploratory"], "^pvalue < 0.05")
})

test_that("run-wide pvalue_column default applies when grid has no p_column", {
  res <- data.frame(
    gene_name = c("A", "B"),
    pvalue = c(0.01, 0.20),
    padj   = c(0.30, 0.90),
    log2FoldChange = c(2, 2),
    stringsAsFactors = FALSE
  )
  grid <- data.frame(name = "standard", p_cutoff = 0.05, log2fc = 1,
                     stringsAsFactors = FALSE)
  sets <- build_deg_threshold_sets(list(T_vs_C = res), grid, pvalue_column = "pvalue")
  expect_equal(sum(sets$standard$T_vs_C$significance != "Not_Sig"), 1)
})
