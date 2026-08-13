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

test_that("report coverage makes missing and covered domains explicit", {
  root <- tempfile()
  dir.create(file.path(root, "1-DEG"), recursive = TRUE)
  dir.create(file.path(root, "3-Visualization"), recursive = TRUE)
  write.csv(data.frame(sample = c("S1", "S2"), condition = c("A", "B")),
            file.path(root, "1-DEG", "colData.csv"), row.names = FALSE)
  file.create(file.path(root, "3-Visualization", "PCA_plot.pdf"))
  file.create(file.path(root, "1-DEG", "DEG_threshold_summary.csv"))
  file.create(file.path(root, "sessionInfo.txt"))
  manifest <- build_report_coverage_manifest(root)
  expect_equal(manifest$status[manifest$domain == "sample_qc"], "covered")
  expect_equal(manifest$report_location[manifest$domain == "sample_qc"],
               "main_body_or_appendix_index")
  expect_equal(manifest$status[manifest$domain == "gsea"], "not_available")
  expect_error(validate_report_coverage(manifest), "Required report domains")
  unlink(root, recursive = TRUE)
})

test_that("claim-evidence ledger requires calibrated levels and provenance", {
  ledger <- data.frame(
    claim_id = "C1", claim = "A differed from B", claim_level = "E1_association",
    scope = "in vitro", evidence = "DE result", source_files = "1-DEG/result.csv",
    assumptions = "valid design", alternatives = "batch", falsifier = "independent null",
    next_evidence = "repeat experiment"
  )
  expect_true(validate_claim_evidence_ledger(ledger))
  ledger$claim_level <- "mechanism_proved"
  expect_error(validate_claim_evidence_ledger(ledger), "unsupported evidence levels")
})
