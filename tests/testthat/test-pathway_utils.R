# Tests for pathway_utils.R helpers
#
# These tests cover both data transformations and focused PDF rendering because
# sample alignment, missing-value colors, and annotation layout can otherwise
# break silently even when the numerical helpers still pass.

# --- internal helpers --------------------------------------------------------

test_that("format_p_stars maps adjusted p-values to star labels", {
  expect_equal(format_p_stars(c(0.0005, 0.005, 0.03, 0.2, NA)),
               c("***", "**", "*", "ns", "ns"))
})

# --- melt_gene_expression ----------------------------------------------------

test_that("melt_gene_expression returns one row per gene x sample", {
  expr <- matrix(1:24, nrow = 4,
                 dimnames = list(paste0("G", 1:4), paste0("S", 1:6)))
  grp <- rep(c("A", "B"), each = 3)
  long <- melt_gene_expression(expr, c("G1", "G3"), grp)
  expect_named(long, c("Gene", "sample", "value", "condition"))
  expect_equal(nrow(long), 2 * 6)
  expect_setequal(unique(long$Gene), c("G1", "G3"))
  # values line up with the source matrix
  expect_equal(long$value[long$Gene == "G1"], as.numeric(expr["G1", ]))
})

test_that("melt_gene_expression drops genes absent from the matrix", {
  expr <- matrix(1:12, nrow = 3, dimnames = list(paste0("G", 1:3), paste0("S", 1:4)))
  long <- melt_gene_expression(expr, c("G1", "NOPE"), rep("A", 4))
  expect_setequal(unique(long$Gene), "G1")
  expect_equal(nrow(melt_gene_expression(expr, "NOPE", rep("A", 4))), 0)
})

test_that("melt_gene_expression aligns a named group vector by sample ID", {
  expr <- matrix(1:8, nrow = 2,
                 dimnames = list(c("G1", "G2"), c("S1", "S2", "S3", "S4")))
  grp <- c(S3 = "B", S1 = "A", S4 = "B", S2 = "A")
  long <- melt_gene_expression(expr, "G1", grp)
  expect_equal(long$condition, c("A", "A", "B", "B"))
  expect_error(melt_gene_expression(expr, "G1", grp[-1]), "must have length")
})

# --- score_gene_sets ---------------------------------------------------------

test_that("score_gene_sets returns a gene_set x sample matrix", {
  testthat::skip_if_not_installed("GSVA")
  expr <- matrix(stats::rnorm(50 * 6), nrow = 50,
                 dimnames = list(paste0("G", 1:50), paste0("S", 1:6)))
  gs <- list(P1 = paste0("G", 1:10), P2 = paste0("G", 11:30))
  sc <- score_gene_sets(expr, gs, method = "zscore")
  expect_equal(dim(sc), c(2, 6))
  expect_equal(rownames(sc), c("P1", "P2"))
  expect_equal(colnames(sc), paste0("S", 1:6))
  audit <- attr(sc, "gene_set_audit")
  expect_equal(audit$overlap_size, c(10L, 20L))
  expect_true(all(audit$retained))
})

test_that("score_gene_sets rejects malformed matrices and undersized overlap", {
  testthat::skip_if_not_installed("GSVA")
  expr <- matrix(stats::rnorm(10 * 4), nrow = 10,
                 dimnames = list(paste0("G", 1:10), paste0("S", 1:4)))
  expect_error(score_gene_sets(expr, list(tiny = c("G1", "G2")), method = "zscore"),
               "No gene set")
  bad <- expr
  rownames(bad)[2] <- rownames(bad)[1]
  expect_error(score_gene_sets(bad, list(ok = paste0("G", 1:6)), method = "zscore"),
               "unique")
})

# --- pathway_group_comparison ------------------------------------------------

test_that("pathway_group_comparison reports delta, direction and global BH", {
  # Two pathways, three groups of n=3 with a clear shift in P1 only.
  set.seed(42)
  grp <- factor(rep(c("PBS", "IZER", "TES"), each = 3), levels = c("PBS", "IZER", "TES"))
  sc <- rbind(
    P1 = c(rnorm(3, 0, 0.1), rnorm(3, 2, 0.1), rnorm(3, 0, 0.1)),
    P2 = c(rnorm(3, 0, 0.1), rnorm(3, 0, 0.1), rnorm(3, 0, 0.1))
  )
  colnames(sc) <- paste0(grp, "_", rep(1:3, 3))
  cmp <- pathway_group_comparison(sc, grp, c("PBS", "IZER", "TES"),
                                  comparisons = list(c("PBS", "IZER"), c("TES", "IZER")))
  expect_named(cmp, c("Pathway", "Comparison", "group1", "group2", "mean1",
                      "mean2", "delta", "p", "p.adj", "direction"))
  # delta = mean(treatment) - mean(control); Comparison = treatment_vs_control
  p1 <- cmp[cmp$Pathway == "P1" & cmp$Comparison == "IZER_vs_PBS", ]
  expect_equal(p1$delta, p1$mean2 - p1$mean1)
  expect_equal(p1$direction, "UP")
  expect_true(p1$p.adj < 0.05)
  # global adjustment: p.adj is a single BH pass over all rows, not per-pathway
  expect_equal(cmp$p.adj, stats::p.adjust(cmp$p, method = "BH"))
})

test_that("pathway_group_comparison returns empty df when no comparisons possible", {
  sc <- matrix(1:6, nrow = 1, dimnames = list("P1", paste0("S", 1:6)))
  grp <- factor(rep("A", 6))
  expect_equal(nrow(pathway_group_comparison(sc, grp, comparisons = list(c("A", "B")))), 0)
})

test_that("pathway_group_comparison reorders named groups and distinguishes zero change", {
  sc <- rbind(P1 = c(0, 0, 1, 1), P2 = c(2, 2, 2, 2))
  colnames(sc) <- c("A1", "A2", "B1", "B2")
  shuffled <- c(B2 = "B", A1 = "A", B1 = "B", A2 = "A")
  cmp <- pathway_group_comparison(sc, shuffled, c("A", "B"),
                                  comparisons = list(c("A", "B")))
  expect_equal(cmp$delta[cmp$Pathway == "P1"], 1)
  expect_equal(cmp$direction[cmp$Pathway == "P2"], "NO_CHANGE")
  expect_error(pathway_group_comparison(sc, c(A1 = "A", A2 = "A", X = "B", B2 = "B"),
                                        c("A", "B"), comparisons = list(c("A", "B"))),
               "does not match")
})

test_that("pathway summary and key-gene heatmap render focused PDFs", {
  testthat::skip_if_not_installed("ComplexHeatmap")
  testthat::skip_if_not_installed("circlize")
  tmp_delta <- tempfile(fileext = ".pdf")
  tmp_heat <- tempfile(fileext = ".pdf")
  on.exit(unlink(c(tmp_delta, tmp_heat)), add = TRUE)

  comp <- data.frame(
    Pathway = c("P1", "P2"), Comparison = c("B_vs_A", "B_vs_A"),
    delta = c(-1, 1), p.adj = c(0.2, 0.01)
  )
  plot_pathway_delta_summary_pdf(comp, tmp_delta)
  expect_true(file.exists(tmp_delta) && file.info(tmp_delta)$size > 0)

  fc <- data.frame(
    Gene = c("G1", "G1", "G2"),
    Comparison = c("B_vs_A", "C_vs_A", "B_vs_A"),
    log2FC = c(1, NA, -1)
  )
  ht <- plot_keygenes_log2fc_heatmap_pdf(
    fc, tmp_heat, row_group = c("Set1", "Set1", "Set2"),
    row_group_levels = c("Set1", "Set2")
  )
  expect_true(file.exists(tmp_heat) && file.info(tmp_heat)$size > 0)
  expect_true(anyNA(ht@matrix))

  dup <- rbind(fc, fc[1, ])
  expect_error(plot_keygenes_log2fc_heatmap_pdf(dup, tempfile(fileext = ".pdf")),
               "duplicate gene-by-comparison")
})
