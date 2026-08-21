test_that("narrative scaffold accepts an explicit template outside project git", {
  project <- tempfile("project-")
  dir.create(project)
  template <- tempfile(fileext = ".Rmd")
  writeLines("# Narrative", template)
  out <- scaffold_narrative_report(project, template = template)
  expect_true(file.exists(out))
  expect_identical(readLines(out), "# Narrative")
  expect_error(scaffold_narrative_report(project, template = template), "Refusing to overwrite")
  unlink(project, recursive = TRUE)
  unlink(template)
})

test_that("make_report_citers reads DEG/gene/GSEA/ORA values from a run layout", {
  run <- tempfile("run-")
  dir.create(file.path(run, "1-DEG", "all_genes"), recursive = TRUE)
  dir.create(file.path(run, "2-GSEA", "standard"), recursive = TRUE)

  write.csv(
    data.frame(Threshold = c("standard", "exploratory"),
               Comparison = c("A_vs_B", "A_vs_B"),
               Total = c(12, 40)),
    file.path(run, "1-DEG", "DEG_threshold_summary.csv"), row.names = FALSE)
  write.csv(
    data.frame(gene_name = c("Gzma", "Cdh1"),
               log2FoldChange_shrunken = c(2.34, -1.21),
               padj = c(0.0004, 0.03)),
    file.path(run, "1-DEG", "all_genes", "DESeq2_all_genes_A_vs_B.csv"), row.names = FALSE)
  write.csv(
    data.frame(Description = c("ribosome biogenesis", "cell killing"),
               NES = c(3.31, 2.05)),
    file.path(run, "2-GSEA", "GSEA_GO_A_vs_B.csv"), row.names = FALSE)
  write.csv(
    data.frame(Description = c("term z", "term a"), p.adjust = c(0.05, 0.001)),
    file.path(run, "2-GSEA", "standard", "GO_ORA_UP_A_vs_B.csv"), row.names = FALSE)

  h <- make_report_citers(run)
  expect_equal(h$sv("standard", "A_vs_B"), "12")
  expect_equal(h$sv("exploratory", "A_vs_B"), "40")
  expect_equal(h$gl("Gzma", "A_vs_B"), "+2.3")
  expect_equal(h$gp("Gzma", "A_vs_B"), "4.00e-04")
  expect_equal(h$gp("Cdh1", "A_vs_B"), "0.030")
  expect_equal(h$gnes("ribosome", "A_vs_B"), "+3.31")
  expect_equal(h$ora_terms("A_vs_B", "UP", "GO", n = 1), "term a")  # ordered by p.adjust

  # Missing values degrade to NA / placeholder rather than erroring.
  expect_equal(h$sv("standard", "MISSING"), "NA")
  expect_true(is.na(h$deg_val("NoGene", "A_vs_B")))
  expect_equal(h$gnes("nomatch", "A_vs_B"), "NA")
  expect_equal(h$ora_terms("A_vs_B", "DOWN", "GO"), "（无显著条目）")
})

test_that("validate_narrative_report passes a good report and fails broken ones", {
  good <- tempfile(fileext = ".html")
  writeLines(c("<h1>技术摘要</h1>", "<h1>差异表达</h1>",
               '<img src="data:image/png;base64,AAAA">',
               '<img src="data:image/png;base64,BBBB">'), good)

  res <- validate_narrative_report(good, min_figures = 2,
                                   required_sections = c("技术摘要", "差异表达"))
  expect_true(res$passed)
  expect_equal(res$figures, 2)

  expect_error(validate_narrative_report(good, min_figures = 5), "embedded figures")
  expect_error(validate_narrative_report(good, required_sections = "综合"), "Required section")
  expect_error(validate_narrative_report(file.path(tempdir(), "nope.html")), "not found")

  missing_fig <- tempfile(fileext = ".html")
  writeLines(c('<img src="data:image/png;base64,AAAA">', "_图不可用：foo.pdf_"), missing_fig)
  expect_error(validate_narrative_report(missing_fig), "Missing-figure marker")

  na_report <- tempfile(fileext = ".html")
  writeLines(c('<img src="data:image/png;base64,AAAA">', "<p>NES ≈ NA</p>"), na_report)
  expect_message(validate_narrative_report(na_report), "NA")
})
