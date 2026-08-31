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

test_that("voom sample fraction applies a total-sample raw-count requirement", {
  skip_if_not_installed("edgeR")
  counts <- rbind(
    ubiquitous = rep(100, 6),
    three_samples = c(30, 30, 30, 0, 0, 0),
    four_samples = c(30, 30, 30, 30, 0, 0),
    low_counts = rep(1, 6)
  )
  colnames(counts) <- paste0("S", seq_len(ncol(counts)))
  group <- rep(c("A", "B"), each = 3)
  half <- prepare_dge_for_voom(counts, group, min_sample_frac = 0.5)
  two_thirds <- prepare_dge_for_voom(counts, group, min_sample_frac = 2 / 3)
  all <- prepare_dge_for_voom(counts, group, min_sample_frac = 1)
  expect_setequal(rownames(half), c("ubiquitous", "three_samples", "four_samples"))
  expect_setequal(rownames(two_thirds), c("ubiquitous", "four_samples"))
  expect_identical(rownames(all), "ubiquitous")
  expect_equal(two_thirds$samples$lib.size, unname(colSums(two_thirds$counts)))
})

test_that("zero voom sample fraction preserves edgeR group-aware filtering", {
  skip_if_not_installed("edgeR")
  counts <- rbind(
    ubiquitous = rep(100, 6),
    small_group = c(30, 30, 0, 0, 0, 0),
    low_counts = rep(1, 6)
  )
  colnames(counts) <- paste0("S", seq_len(ncol(counts)))
  group <- c("A", "A", rep("B", 4))
  dge <- edgeR::DGEList(counts, group = group)
  expected <- edgeR::filterByExpr(dge, min.count = 10)
  legacy <- prepare_dge_for_voom(counts, group, min_sample_frac = 0)
  half <- prepare_dge_for_voom(counts, group, min_sample_frac = 0.5)
  expect_identical(rownames(legacy), rownames(counts)[expected])
  expect_true("small_group" %in% rownames(legacy))
  expect_false("small_group" %in% rownames(half))
})

test_that("meeting the raw-count fraction cannot bypass edgeR library-size filtering", {
  skip_if_not_installed("edgeR")
  counts <- rbind(
    ubiquitous = rep(100, 6),
    large_libraries = c(1000000, 1000000, 0, 0, 0, 0),
    raw_fraction_only = c(10, 10, 10, 0, 0, 0)
  )
  colnames(counts) <- paste0("S", seq_len(ncol(counts)))
  group <- rep(c("A", "B"), each = 3)
  expect_equal(sum(counts["raw_fraction_only", ] >= 10), 3L)
  filtered <- prepare_dge_for_voom(counts, group, min_sample_frac = 0.5)
  expect_identical(rownames(filtered), "ubiquitous")
})

test_that("voom filter parameters are validated and empty filtering is explicit", {
  skip_if_not_installed("edgeR")
  counts <- matrix(20, 4, 6)
  group <- rep(c("A", "B"), each = 3)
  for (fraction in list(-0.1, 1.1, NA_real_, Inf, c(0.3, 0.5))) {
    expect_error(prepare_dge_for_voom(counts, group, min_sample_frac = fraction),
                 "min_sample_frac")
  }
  expect_error(prepare_dge_for_voom(counts, group, min_counts_per_sample = -1),
               "min_counts_per_sample")
  expect_error(prepare_dge_for_voom(counts, group, min_counts_per_sample = 100),
               "No genes remain")
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
