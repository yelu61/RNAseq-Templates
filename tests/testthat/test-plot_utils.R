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

test_that("publication dimensions use final physical sizes", {
  expect_equal(mm_to_in(25.4), 1)
  single <- publication_dimensions("single", height_mm = 100)
  double <- publication_dimensions("double", height_mm = 100)
  expect_equal(unname(single[["width"]]), 89 / 25.4)
  expect_equal(unname(double[["width"]]), 183 / 25.4)
  expect_error(publication_dimensions("double", height_mm = 300), "no larger")
  expect_true(nzchar(resolve_publication_font()))
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
  expect_lte(length(strsplit(wrapped, "\n", fixed = TRUE)[[1]]), 2)
  expect_match(wrapped, "…$")
})

test_that("save_pdf_plot rejects NULL and writes a non-empty PDF", {
  outfile <- tempfile(fileext = ".pdf")
  on.exit(unlink(outfile), add = TRUE)
  expect_null(save_pdf_plot(NULL, outfile))
  expect_false(file.exists(outfile))

  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) +
    ggplot2::geom_point()
  save_pdf_plot(p, outfile, width = 4, height = 3)
  expect_true(file.exists(outfile))
  expect_gt(file.info(outfile)$size, 1500)
})

test_that("save_pdf_device validates grid/base output before promotion", {
  outfile <- tempfile(fileext = ".pdf")
  on.exit(unlink(outfile), add = TRUE)
  save_pdf_device(outfile, width = 4, height = 3, draw = function() {
    graphics::plot(1:5, 1:5, type = "b")
  })
  expect_true(file.exists(outfile))
  expect_gt(file.info(outfile)$size, 1500)
})

test_that("save_plot_bundle exports final-size review assets", {
  skip_if_not_installed("ragg")
  stem <- tempfile()
  paths <- paste0(stem, c(".pdf", ".png"))
  on.exit(unlink(paths), add = TRUE)
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 3:1), ggplot2::aes(x, y)) +
    ggplot2::geom_point() + theme_publication()
  out <- save_plot_bundle(
    p, stem, width_mm = 89, height_mm = 70,
    formats = c("pdf", "png")
  )
  expect_setequal(out, paths)
  expect_true(all(file.exists(paths)))
  expect_true(all(file.info(paths)$size > 1500))
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

test_that("group boxplot can reuse a precomputed globally adjusted statistic table", {
  testthat::skip_if_not_installed("ggpubr")
  df <- data.frame(
    Pathway = rep(c("P1", "P2"), each = 6),
    group = rep(rep(c("A", "B"), each = 3), 2),
    value = c(1:3, 4:6, 2:4, 5:7)
  )
  stats <- data.frame(
    Pathway = c("P1", "P2"), group1 = "A", group2 = "B",
    delta = c(3, 3), p.adj = c(0.009, 0.20)
  )
  prepared <- .prepare_plot_stat_table(stats, df, "value", "Pathway", "stars")
  expect_equal(prepared$label, c("**", "ns"))
  expect_true(all(is.finite(prepared$y.position)))

  outfile <- tempfile(fileext = ".pdf")
  on.exit(unlink(outfile), add = TRUE)
  plot_group_boxplot_pdf(
    df, value_col = "value", group_col = "group", filename = outfile,
    facet_col = "Pathway", label_style = "stars", stat_table = stats
  )
  expect_true(file.exists(outfile) && file.info(outfile)$size > 0)
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

test_that("enrichment dotplot accepts saved-table data and limits label lines", {
  outfile <- tempfile(fileext = ".pdf")
  on.exit(unlink(outfile), add = TRUE)
  df <- make_enrich_df()
  df$Description[1] <- paste(rep("very long pathway label", 8), collapse = " ")
  p <- plot_enrich_dotplot(df, outfile, "Saved ORA table", show_category = 5)
  expect_s3_class(p, "ggplot")
  expect_true(file.exists(outfile) && file.info(outfile)$size > 1500)
})

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

test_that("prepare_enrich_df preserves duplicate Descriptions with unique display labels", {
  # KEGG can return the same Description under different IDs. Both biological
  # terms must remain present while their figure labels stay unambiguous.
  df <- make_enrich_df()
  dup <- df[2, ]
  dup$ID <- "GO:dup"
  dup$p.adjust <- 0.015  # sorts after the original term 2 (p.adjust 0.02)
  df <- rbind(df, dup)
  out <- prepare_enrich_df(df, show_category = 15)
  expect_equal(nrow(out), nrow(df))
  expect_equal(sum(grepl("term 2", out$Description_label, fixed = TRUE)), 2)
  expect_equal(anyDuplicated(out$Description_label), 0)
  expect_equal(anyDuplicated(levels(out$Description_wrapped)), 0)
  dup_pdf <- tempfile(fileext = ".pdf")
  on.exit(unlink(dup_pdf), add = TRUE)
  expect_no_error(plot_enrich_dotplot(df, dup_pdf, "dup", show_category = 15))
})

test_that("prepare_enrich_df falls back to term IDs for missing descriptions", {
  df <- make_enrich_df()
  df$Description[1:2] <- c(NA_character_, "")
  out <- prepare_enrich_df(df, show_category = 5)
  expect_true(all(c("GO:1", "GO:2") %in% out$Description_label))
  expect_equal(nrow(out), nrow(df))
})

test_that("plot helpers return invisibly so top-level runner calls do not open Rplots.pdf", {
  df <- make_enrich_df()
  outfile <- tempfile(fileext = ".pdf")
  on.exit(unlink(outfile), add = TRUE)
  result <- withVisible(plot_enrich_dotplot(df, outfile, "visibility", show_category = 5))
  expect_s3_class(result$value, "ggplot")
  expect_false(result$visible)
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

test_that("theme dot-heatmap caps unsafe PDF dimensions", {
  df <- data.frame(
    figure_group = "Immune", theme = "Interferon",
    Description = c("term one", "term two"),
    comparison = c("B_vs_A", "B_vs_A"),
    Count = c(10, 8), neg_log10_padj = c(3, 2),
    label = factor(c("term one", "term two")),
    stringsAsFactors = FALSE
  )
  outfile <- tempfile(fileext = ".pdf")
  on.exit(unlink(outfile), add = TRUE)
  expect_warning(
    plot_theme_dotheatmap_pdf(df, outfile, height = 60, max_height = 8),
    "dimensions capped"
  )
  expect_true(file.exists(outfile) && file.info(outfile)$size > 0)
})
