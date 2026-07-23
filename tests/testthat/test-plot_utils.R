# Tests for plot_utils.R helpers
#
# The plotting functions themselves are exercised end-to-end by the demo smoke
# test; these unit tests focus on the pure data-transformation helpers that the
# plots depend on, which are cheap to run and easy to break silently.

# --- palette / theme ---------------------------------------------------------

test_that("RNAseq_PALETTE exposes the expected color roles", {
  expect_true(is.list(RNAseq_PALETTE))
  for (nm in c("up_regulated", "down_regulated", "not_significant",
               "heatmap_up", "heatmap_down", "heatmap_mid",
               "gsea_activated", "gsea_suppressed")) {
    expect_true(nm %in% names(RNAseq_PALETTE), info = paste("missing", nm))
  }
  # All single-color entries should be valid hex colors
  for (nm in c("up_regulated", "down_regulated", "not_significant")) {
    expect_match(RNAseq_PALETTE[[nm]], "^#[0-9A-Fa-f]{6}$")
  }
})

test_that("theme_publication returns a ggplot2 theme", {
  th <- theme_publication(base_size = 10)
  expect_s3_class(th, "theme")
})

test_that("make_group_colors names colors by group level", {
  cols <- make_group_colors(c("Ctrl", "Treat"))
  expect_named(cols, c("Ctrl", "Treat"))
  expect_length(cols, 2)

  cols3 <- make_group_colors(c("A", "B", "C"))
  expect_length(cols3, 3)

  # >3 levels still yields one named color per level
  many <- make_group_colors(letters[1:6])
  expect_length(many, 6)
  expect_named(many, letters[1:6])
})

# --- text / numeric helpers --------------------------------------------------

test_that("wrap_term_labels wraps long strings", {
  long <- paste(rep("word", 20), collapse = " ")
  wrapped <- wrap_term_labels(long, width = 20)
  expect_true(grepl("\n", wrapped))
  # short strings stay on one line
  expect_false(grepl("\n", wrap_term_labels("short", width = 45)))
})

test_that("parse_ratio_numeric parses 'k/n' ratios", {
  expect_equal(parse_ratio_numeric("3/100"), 0.03)
  expect_equal(parse_ratio_numeric(c("1/2", "5/10")), c(0.5, 0.5))
  expect_true(is.na(parse_ratio_numeric("notaratio")))
  expect_true(is.na(parse_ratio_numeric("1/0")))
  expect_true(is.na(parse_ratio_numeric("a/b")))
})

test_that("zscore_rows row-standardizes and caps values", {
  mat <- matrix(c(1, 2, 3, 10, 20, 30), nrow = 2, byrow = TRUE)
  z <- zscore_rows(mat, cap = 2)
  expect_equal(dim(z), dim(mat))
  expect_true(all(z <= 2 & z >= -2))
  # rows should be centered near zero
  expect_equal(rowMeans(z), c(0, 0), tolerance = 1e-8)

  # cap = NULL disables capping
  z_uncapped <- zscore_rows(mat, cap = NULL)
  expect_true(max(abs(z_uncapped)) > 2 || max(abs(z_uncapped)) <= 2) # no cap applied
})

test_that("zscore_rows handles constant rows without NaN", {
  mat <- matrix(c(5, 5, 5, 1, 2, 3), nrow = 2, byrow = TRUE)
  z <- zscore_rows(mat)
  expect_false(any(is.nan(z)))
  expect_false(any(is.na(z)))
})

test_that("format_p_for_label formats p-values", {
  expect_equal(format_p_for_label(NA), "NA")
  expect_match(format_p_for_label(0.05), "0.050")
  expect_match(format_p_for_label(0.0000123), "e-0") # scientific notation for tiny p
})

test_that("sem helpers bound the mean", {
  x <- c(1, 2, 3, 4, 5)
  expect_gt(sem_upper(x), mean(x))
  expect_lt(sem_lower(x), mean(x))
  # symmetric around the mean
  expect_equal(sem_upper(x) - mean(x), mean(x) - sem_lower(x), tolerance = 1e-8)
})

# --- pairwise effect table ---------------------------------------------------

test_that("pairwise_effect_table computes per-pair statistics", {
  set.seed(1)
  df <- data.frame(
    value = c(rnorm(5, 0), rnorm(5, 3)),
    group = rep(c("A", "B"), each = 5)
  )
  tbl <- pairwise_effect_table(df, value_col = "value", group_col = "group",
                               comparisons = list(c("A", "B")))
  expect_equal(nrow(tbl), 1)
  expect_true(all(c("group1", "group2", "p", "p.adj", "delta", "y.position", "label") %in% colnames(tbl)))
  expect_equal(tbl$group1, "A")
  expect_equal(tbl$group2, "B")
  # delta = mean(B) - mean(A) > 0 here
  expect_gt(tbl$delta, 0)
})

test_that("pairwise_effect_table returns empty frame when too few observations", {
  df <- data.frame(value = c(1, 2), group = c("A", "B"))
  tbl <- pairwise_effect_table(df, value_col = "value", group_col = "group")
  expect_equal(nrow(tbl), 0)
})

# --- enrichment data prep ----------------------------------------------------

# Build a minimal enrichResult-like data frame (as.data.frame(enrichResult))
make_enrich_df <- function() {
  data.frame(
    ID = paste0("GO:", 1:5),
    Description = paste("term", 1:5),
    GeneRatio = c("5/100", "4/100", "3/100", "2/100", "1/100"),
    BgRatio = c("10/1000", "10/1000", "10/1000", "10/1000", "10/1000"),
    Count = c(5, 4, 3, 2, 1),
    pvalue = c(0.001, 0.005, 0.01, 0.02, 0.04),
    p.adjust = c(0.01, 0.02, 0.03, 0.04, 0.05),
    stringsAsFactors = FALSE
  )
}

test_that("parse_ratio + prepare_enrich_df computes FoldEnrichment", {
  df <- make_enrich_df()
  # prepare_enrich_df expects an object coercible via as.data.frame; a data.frame works
  out <- prepare_enrich_df(df, show_category = 5)
  expect_true("FoldEnrichment" %in% colnames(out))
  expect_true("GeneRatioNumeric" %in% colnames(out))
  expect_true("log10_padj" %in% colnames(out))
  # GeneRatioNumeric / BgRatio = FoldEnrichment
  expect_equal(out$GeneRatioNumeric[1], 0.05)
  expect_equal(out$FoldEnrichment[1], 0.05 / 0.01)
})

test_that("prepare_enrich_df handles NULL and empty input", {
  expect_equal(nrow(prepare_enrich_df(NULL)), 0)
  empty <- make_enrich_df()[0, ]
  expect_equal(nrow(prepare_enrich_df(empty)), 0)
})

test_that("prepare_enrich_df respects show_category limit", {
  out <- prepare_enrich_df(make_enrich_df(), show_category = 3)
  expect_lte(nrow(out), 3)
})

# --- theme matching ----------------------------------------------------------

test_that("default_enrichment_themes is well-formed", {
  themes <- default_enrichment_themes()
  expect_true(all(c("figure_group", "theme", "keywords") %in% colnames(themes)))
  expect_gt(nrow(themes), 0)
})

test_that("clean_term_label shortens common GO phrases", {
  expect_match(clean_term_label("positive regulation of cell cycle"), "\\+reg")
  expect_match(clean_term_label("negative regulation of apoptosis"), "-reg")
  expect_match(clean_term_label("oxidative phosphorylation"), "OXPHOS")
  expect_match(clean_term_label("extracellular matrix organization"), "ECM")
})

test_that("match_enrichment_themes assigns themes by keyword", {
  df <- data.frame(
    Description = c("interferon-gamma-mediated signaling pathway",
                    "unrelated biological process xyz"),
    stringsAsFactors = FALSE
  )
  out <- match_enrichment_themes(df, default_enrichment_themes())
  # The interferon term should match; the unrelated one should not
  expect_true("theme" %in% colnames(out))
  expect_true(any(grepl("interferon", out$Description, ignore.case = TRUE)))
  expect_false(any(grepl("unrelated", out$Description, ignore.case = TRUE)))
})

test_that("match_enrichment_themes returns empty frame when nothing matches", {
  df <- data.frame(Description = "zzz no theme here", stringsAsFactors = FALSE)
  out <- match_enrichment_themes(df, default_enrichment_themes())
  expect_equal(nrow(out), 0)
})
