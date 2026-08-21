test_that("run-review notebook is read-only with respect to core analysis", {
  root <- if (exists("repo_root", inherits = TRUE)) {
    get("repo_root", inherits = TRUE)
  } else {
    rprojroot::find_root(rprojroot::is_git_root, path = getwd())
  }
  path <- file.path(root, "notebooks", "RNAseq_General_RunReview.R")
  expect_true(file.exists(path))
  lines <- readLines(path, warn = FALSE)
  code <- lines[!grepl("^#", lines)]

  forbidden <- c("DESeqDataSet", "DESeq\\(", "lfcShrink\\(",
                 "enrichGO\\(", "enrichKEGG\\(", "gseGO\\(", "gseKEGG\\(")
  for (pattern in forbidden) expect_false(any(grepl(pattern, code)))
  expect_true(any(grepl("REVIEW_OUTDIR must be outside RUN_DIR", lines, fixed = TRUE)))
  expect_true(any(grepl("review_provenance.csv", lines, fixed = TRUE)))
})
