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

test_that("run_voom contains plotting inside the requested PDF device", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("limma")
  counts <- matrix(
    c(50, 55, 52, 100, 110, 105,
      80, 75, 78, 40, 45, 42,
      30, 32, 31, 60, 63, 61,
      90, 88, 92, 95, 96, 94),
    nrow = 4, byrow = TRUE,
    dimnames = list(paste0("g", 1:4), paste0("s", 1:6))
  )
  group <- factor(rep(c("A", "B"), each = 3))
  dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts, group = group))
  design <- make_group_design(group)
  outfile <- tempfile(fileext = ".pdf")
  on.exit(unlink(outfile), add = TRUE)
  old_device <- getOption("device")
  on.exit(options(device = old_device), add = TRUE)
  options(device = function(...) stop("unexpected default graphics device"))
  v <- run_voom(dge, design = design, plot_file = outfile)
  expect_true(is.matrix(v$E))
  expect_true(file.exists(outfile) && file.info(outfile)$size > 1500)
})
