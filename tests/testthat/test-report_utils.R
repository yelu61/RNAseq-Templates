# Tests for report_utils.R (figure listing and safe CSV reading)

test_that("read_csv_safe returns NULL for missing file and a df for real file", {
  expect_null(read_csv_safe("/nonexistent/path/file.csv"))

  tmp <- tempfile(fileext = ".csv")
  write.csv(data.frame(a = 1:3, b = letters[1:3]), tmp, row.names = FALSE)
  out <- read_csv_safe(tmp)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 3)
  unlink(tmp)
})

test_that("read_csv_safe returns NULL on parse error instead of throwing", {
  tmp <- tempfile(fileext = ".csv")
  # A directory with a .csv name causes read.csv to error.
  dir.create(tmp)
  # read.csv emits warnings before erroring; suppress them — we only care that
  # read_csv_safe swallows the failure and returns NULL rather than throwing.
  expect_null(suppressWarnings(read_csv_safe(tmp)))
  unlink(tmp, recursive = TRUE)
})

test_that("list_report_figures finds PDFs recursively and returns relative paths", {
  base <- tempfile()
  dir.create(file.path(base, "3-Visualization", "sub"), recursive = TRUE)
  pdf1 <- file.path(base, "3-Visualization", "a.pdf")
  pdf2 <- file.path(base, "3-Visualization", "sub", "b.pdf")
  file.create(pdf1, pdf2)
  file.create(file.path(base, "3-Visualization", "ignore.txt"))

  figs <- list_report_figures(file.path(base, "3-Visualization"), base = base)
  expect_length(figs, 2)
  expect_true(all(grepl("\\.pdf$", figs)))
  expect_false(any(grepl("ignore", figs)))
  # Relative to base, not absolute.
  expect_false(any(grepl(base, figs, fixed = TRUE)))
  unlink(base, recursive = TRUE)
})

test_that("list_report_figures returns empty for missing dir", {
  expect_equal(list_report_figures("/nonexistent/dir"), character(0))
})
