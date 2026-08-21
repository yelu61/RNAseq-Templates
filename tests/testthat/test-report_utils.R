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

test_that("coverage overrides cannot replace refreshed evidence inventory", {
  root <- tempfile()
  dir.create(file.path(root, "2-GSEA"), recursive = TRUE)
  writeLines("term", file.path(root, "2-GSEA", "GSEA_GO_case.csv"))
  overrides <- data.frame(
    domain = "gsea", status = "omitted_with_reason",
    report_location = "appendix", omission_reason = "Not in headline report",
    evidence_count = 0, evidence_files = "stale.csv", stringsAsFactors = FALSE
  )
  manifest <- build_report_coverage_manifest(root, overrides = overrides)
  row <- manifest[manifest$domain == "gsea", ]
  expect_equal(row$status, "omitted_with_reason")
  expect_equal(row$evidence_count, 1)
  expect_match(row$evidence_files, "GSEA_GO_case.csv", fixed = TRUE)
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

test_that("claim-evidence ledger verifies project-relative source files", {
  root <- tempfile()
  dir.create(file.path(root, "1-DEG"), recursive = TRUE)
  writeLines("evidence", file.path(root, "1-DEG", "result.csv"))
  ledger <- data.frame(
    claim_id = "C1", claim = "A differed from B", claim_level = "E1_association",
    scope = "in vitro", evidence = "DE result", source_files = "1-DEG/result.csv",
    assumptions = "valid design", alternatives = "batch", falsifier = "independent null",
    next_evidence = "repeat experiment"
  )
  expect_true(validate_claim_evidence_ledger(ledger, base_dir = root))
  ledger$source_files <- "1-DEG/missing.csv"
  expect_error(validate_claim_evidence_ledger(ledger, base_dir = root), "missing source files")
  unlink(root, recursive = TRUE)
})

test_that("interpretation scaffold creates per-section files without overwriting", {
  root <- tempfile()
  dir.create(root)

  res <- scaffold_report_interpretation(root)
  sections <- report_interpretation_sections()
  expect_true(dir.exists(file.path(root, "report_interpretation")))
  expect_identical(sort(res$created), sort(names(sections)))
  for (key in names(sections)) {
    expect_true(file.exists(file.path(root, "report_interpretation", paste0(key, ".md"))))
  }
  # Starter files contain only whole-line HTML comments -> render nothing.
  starter <- readLines(file.path(root, "report_interpretation", "deg.md"), warn = FALSE)
  expect_false(any(nzchar(trimws(strip_interpretation_comments(starter)))))

  # A second run skips everything and preserves written annotations.
  deg_path <- file.path(root, "report_interpretation", "deg.md")
  writeLines("Treatment induces interferon signalling.", deg_path)
  res2 <- scaffold_report_interpretation(root)
  expect_length(res2$created, 0)
  expect_identical(sort(res2$skipped), sort(names(sections)))
  expect_identical(readLines(deg_path, warn = FALSE), "Treatment induces interferon signalling.")

  # overwrite = TRUE regenerates starters.
  res3 <- scaffold_report_interpretation(root, overwrite = TRUE)
  expect_identical(sort(res3$created), sort(names(sections)))
  expect_false(any(grepl("interferon", readLines(deg_path, warn = FALSE))))
  unlink(root, recursive = TRUE)
})

test_that("strip_interpretation_comments removes only whole-line HTML comments", {
  lines <- c(
    "<!-- guidance -->",
    "  <!-- indented guidance -->",
    "Real prose stays.",
    "Inline <!-- comment --> stays too.",
    "<!-- unterminated never matches"
  )
  out <- strip_interpretation_comments(lines)
  expect_identical(out, lines[3:5])
})

test_that("report review scaffold is durable and publish sign-off is strict", {
  root <- tempfile()
  dir.create(root)
  path <- scaffold_report_review(root)
  review <- read.csv(path, stringsAsFactors = FALSE)
  expect_true(validate_report_review(review, require_signoff = FALSE))
  expect_error(validate_report_review(review, require_signoff = TRUE), "remains pending")

  review$status <- "approved"
  review$reviewer <- "reviewer"
  review$reviewed_at <- "2026-08-21"
  write.csv(review, path, row.names = FALSE)
  expect_true(validate_report_review(path, require_signoff = TRUE))

  original <- readLines(path, warn = FALSE)
  scaffold_report_review(root)
  expect_identical(readLines(path, warn = FALSE), original)
  unlink(root, recursive = TRUE)
})

test_that("analysis report validation separates draft rendering from publication", {
  root <- tempfile()
  dir.create(file.path(root, "0-Config"), recursive = TRUE)
  writeLines("config", file.path(root, "0-Config", "analysis_config_used.R"))
  writeLines("session", file.path(root, "sessionInfo.txt"))
  writeLines(paste(rep("x", 11000), collapse = ""), file.path(root, "RNAseq_report.html"))
  dir.create(file.path(root, "1-DEG"), recursive = TRUE, showWarnings = FALSE)
  write.csv(data.frame(Threshold = "standard", Comparison = "A_vs_B", UP = 2,
                       DOWN = 1, Total = 3),
            file.path(root, "1-DEG", "DEG_threshold_summary.csv"), row.names = FALSE)
  writeLines("   - A_vs_B: 3 DEGs (Up: 2, Down: 1)", file.path(root, "Analysis_summary.txt"))

  coverage <- build_report_coverage_manifest(root)
  coverage$status <- "covered"
  coverage$report_location <- "main_body_or_appendix_index"
  coverage$omission_reason <- ""
  write.csv(coverage, file.path(root, "report_coverage_manifest.csv"), row.names = FALSE)
  checklist <- scaffold_report_review(root)

  expect_true(all(validate_analysis_report(root)$passed))
  expect_error(validate_analysis_report(root, publish = TRUE), "not ready to publish")

  review <- read.csv(checklist, stringsAsFactors = FALSE)
  review$status <- "approved"
  review$reviewer <- "reviewer"
  review$reviewed_at <- "2026-08-21"
  write.csv(review, checklist, row.names = FALSE)
  expect_true(all(validate_analysis_report(root, publish = TRUE)$passed))
  unlink(root, recursive = TRUE)
})

test_that("headline validation catches a stale text summary", {
  root <- tempfile()
  dir.create(file.path(root, "1-DEG"), recursive = TRUE)
  write.csv(data.frame(Threshold = "standard", Comparison = "A_vs_B", UP = 2,
                       DOWN = 1, Total = 3),
            file.path(root, "1-DEG", "DEG_threshold_summary.csv"), row.names = FALSE)
  writeLines("   - A_vs_B: 4 DEGs (Up: 3, Down: 1)", file.path(root, "Analysis_summary.txt"))
  expect_error(validate_report_headline_values(root), "disagree")
  unlink(root, recursive = TRUE)
})

test_that("HTML asset validation catches failed previews", {
  html <- tempfile(fileext = ".html")
  writeLines("<html>PDF preview failed: x.pdf</html>", html)
  expect_error(validate_report_html_assets(html), "failure markers")
  writeLines('<html><img src="data:image/png;base64,AAAA"></html>', html)
  out <- validate_report_html_assets(html, expect_figures = TRUE)
  expect_equal(out$embedded_figures, 1L)
  unlink(html)
})

test_that("rmarkdown twin of the .qmd report template exists and does not drift", {
  repo <- tryCatch(rprojroot::find_root(rprojroot::is_git_root), error = function(e) getwd())
  qmd <- file.path(repo, "reports", "analysis_report.qmd")
  rmd <- file.path(repo, "reports", "analysis_report.Rmd")
  expect_true(file.exists(qmd))
  expect_true(file.exists(rmd))
  strip_yaml <- function(path) {
    lines <- readLines(path, warn = FALSE)
    stopifnot(lines[1] == "---")
    lines[-seq_len(which(lines == "---")[2])]
  }
  # Bodies must stay identical; only the YAML front matter differs by design.
  expect_identical(strip_yaml(rmd), strip_yaml(qmd))
})

test_that("render_analysis_report falls back to the .Rmd twin without Quarto", {
  skip_if(requireNamespace("quarto", quietly = TRUE) && nzchar(Sys.which("quarto")),
          "Quarto available: fallback path not exercised")
  skip_if_not(requireNamespace("rmarkdown", quietly = TRUE), "rmarkdown not installed")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")
  outdir <- tempfile()
  dir.create(outdir)
  for (d in c("1-DEG", "2-GSEA", "3-Visualization")) dir.create(file.path(outdir, d))
  report <- file.path(outdir, "RNAseq_report.html")
  expect_message(
    render_analysis_report(outdir = outdir, report_file = report),
    "rmarkdown twin"
  )
  expect_true(file.exists(report))
  expect_gt(file.info(report)$size, 10000)
  unlink(outdir, recursive = TRUE)
})

test_that("report honors the declared primary threshold for detailed DEG tables", {
  skip_if(requireNamespace("quarto", quietly = TRUE) && nzchar(Sys.which("quarto")),
          "Quarto available: fallback path not exercised")
  skip_if_not(requireNamespace("rmarkdown", quietly = TRUE), "rmarkdown not installed")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")
  outdir <- tempfile()
  dir.create(outdir)
  for (d in c("1-DEG/strict", "1-DEG/standard", "2-GSEA", "3-Visualization")) {
    dir.create(file.path(outdir, d), recursive = TRUE)
  }
  deg <- data.frame(
    gene_name = c("StrictGene", "OtherGene"),
    log2FoldChange_shrunken = c(2, -1.5),
    padj = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )
  utils::write.csv(deg, file.path(outdir, "1-DEG/strict/DEG_results_case.csv"), row.names = FALSE)
  deg$gene_name <- c("StandardGene", "OtherGene")
  utils::write.csv(deg, file.path(outdir, "1-DEG/standard/DEG_results_case.csv"), row.names = FALSE)
  report <- file.path(outdir, "RNAseq_report.html")
  render_analysis_report(
    outdir = outdir, report_file = report,
    params = list(primary_threshold = "strict")
  )
  html <- paste(readLines(report, warn = FALSE), collapse = "\n")
  expect_match(html, "StrictGene")
  expect_false(grepl("StandardGene", html, fixed = TRUE))
  unlink(outdir, recursive = TRUE)
})
