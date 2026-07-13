test_that("make_group_design includes batch in the fitted model", {
  group <- factor(c("A", "A", "B", "B"), levels = c("A", "B"))
  batch <- c("X", "Y", "X", "Y")
  design <- make_group_design(group, batch = batch)
  expect_true(all(c("A", "B") %in% colnames(design)))
  expect_true(any(grepl("batch", colnames(design))))
  expect_equal(qr(design)$rank, ncol(design))
})

test_that("make_group_design rejects confounded group and batch", {
  group <- factor(c("A", "A", "B", "B"), levels = c("A", "B"))
  expect_error(make_group_design(group, batch = c("X", "X", "Y", "Y")), "not full rank")
})

test_that("make_group_design rejects invalid batch vectors", {
  group <- factor(c("A", "A", "B", "B"), levels = c("A", "B"))
  expect_error(make_group_design(group, batch = c("X", "Y")), "length")
  expect_error(make_group_design(group, batch = rep("X", 4)), "at least two levels")
})
