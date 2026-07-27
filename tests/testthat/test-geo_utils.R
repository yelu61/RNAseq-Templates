# Tests for geo_utils.R (offline parsing/collapsing only; network download not tested)

test_that("prepare_geo_counts collapses probes to symbols by max mean expression", {
  skip_if_not_installed("dplyr")
  expr <- data.frame(
    S1 = c(1, 5, 2, 9),
    S2 = c(2, 6, 3, 8),
    row.names = c("p1", "p2", "p3", "p4"),
    check.names = FALSE
  )
  feature <- data.frame(
    ID = c("p1", "p2", "p3", "p4"),
    GeneSymbol = c("GeneA", "GeneA", "GeneB", "GeneB"),
    stringsAsFactors = FALSE
  )
  out <- prepare_geo_counts(expr, feature, probe_col = "ID", symbol_col = "GeneSymbol")

  # Result is ordered by descending mean expression: GeneB (max mean 8.5) first.
  expect_equal(sort(rownames(out)), c("GeneA", "GeneB"))
  # GeneA keeps p2 (mean (5+6)/2=5.5) over p1 ((1+2)/2=1.5); GeneB keeps p4 (8.5) over p3 (2.5)
  expect_equal(unname(out["GeneA", "S1"]), 5)
  expect_equal(unname(out["GeneA", "S2"]), 6)
  expect_equal(unname(out["GeneB", "S2"]), 8)
})

test_that("prepare_geo_counts errors clearly when symbol column missing", {
  skip_if_not_installed("dplyr")
  expr <- data.frame(S1 = 1:2, row.names = c("p1", "p2"))
  feature <- data.frame(ID = c("p1", "p2"), foo = c("x", "y"))
  expect_error(prepare_geo_counts(expr, feature), "symbol column not found|Gene symbol")
})
