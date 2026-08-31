# Tests for enrichment_utils.R helpers

test_that("enrich_result_to_df returns empty df for NULL", {
  out <- enrich_result_to_df(NULL)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0)
})

test_that("enrich_result_to_df returns empty df for empty object", {
  out <- enrich_result_to_df(data.frame())
  expect_equal(nrow(out), 0)
})

test_that("build_multi_comparison_enrich_df handles empty list", {
  out <- build_multi_comparison_enrich_df(list())
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0)
})

test_that("build_multi_comparison_enrich_df binds comparison column", {
  df1 <- data.frame(
    ID = "GO:001",
    Description = "term1",
    p.adjust = 0.01,
    Count = 5,
    stringsAsFactors = FALSE
  )
  df2 <- data.frame(
    ID = "GO:002",
    Description = "term2",
    p.adjust = 0.02,
    Count = 3,
    stringsAsFactors = FALSE
  )
  out <- build_multi_comparison_enrich_df(list(A = df1, B = df2))
  expect_equal(nrow(out), 2)
  expect_equal(out$comparison, c("A", "B"))
})

test_that("build_multi_comparison_enrich_df applies ontology filter", {
  df1 <- data.frame(
    ID = "GO:001",
    Description = "term1",
    p.adjust = 0.01,
    Count = 5,
    ONTOLOGY = "BP",
    stringsAsFactors = FALSE
  )
  df2 <- data.frame(
    ID = "GO:002",
    Description = "term2",
    p.adjust = 0.02,
    Count = 3,
    ONTOLOGY = "MF",
    stringsAsFactors = FALSE
  )
  out <- build_multi_comparison_enrich_df(list(A = df1, B = df2), ontology_filter = "BP")
  expect_equal(nrow(out), 1)
  expect_equal(out$ID, "GO:001")
})

test_that("run_go_gsea is reproducible with a fixed seed", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("AnnotationDbi")

  # Fixed ranked ENTREZ list; the function's internal seed governs permutation.
  set.seed(42)
  ids <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "ENTREZID")
  ids <- sample(ids, 2000)
  ranked <- sort(stats::setNames(stats::rnorm(length(ids)), ids), decreasing = TRUE)

  r1 <- run_go_gsea(ranked, org_db = org.Hs.eg.db::org.Hs.eg.db, ont = "BP", p_cutoff = 0.05)
  r2 <- run_go_gsea(ranked, org_db = org.Hs.eg.db::org.Hs.eg.db, ont = "BP", p_cutoff = 0.05)
  expect_s4_class(r1, "gseaResult")
  expect_s4_class(r2, "gseaResult")
  expect_identical(as.data.frame(r1), as.data.frame(r2))
  df <- as.data.frame(r1)
  expect_true(all(!is.na(df$ID)))
  expect_true(all(df$input_overlap >= 10 & df$input_overlap <= 500))
  expect_equal(df$p.adjust, p.adjust(df$pvalue, "BH"))
})

test_that("run_go_gsea returns NULL for a too-short ranked list", {
  skip_if_not_installed("clusterProfiler")
  out <- run_go_gsea(c("1" = 1, "2" = 2), org_db = NULL, min_size = 10)
  expect_null(out)
})

test_that("GSEA audit retains all-NA and partial-NA rows without inventing IDs", {
  df <- data.frame(ID = c(NA, "GO:partial", "GO:finite", "GO:inf", ""),
    Description = c(NA, "partial", NA, "infinite", "unknown"),
    NES = c(NA, 2, -1, Inf, 3), pvalue = c(NA, NA, 0.01, 0.01, 0.001),
    p.adjust = c(NA, NA, 0.2, 0.02, 0.005))
  audit <- audit_gsea_table(df)
  expect_equal(nrow(audit$table), nrow(df))
  expect_identical(audit$table$ID, df$ID)
  expect_equal(audit$summary$valid_terms, 1)
  expect_equal(audit$summary$unusable_rows, 4)
  expect_equal(audit$summary$significant_terms, 0)
  expect_match(audit$table$result_issue[1], "missing_ID")
  expect_match(audit$table$result_issue[2], "nonfinite_pvalue")
  expect_equal(nrow(significant_gsea_terms(df)), 0)
})

test_that("GSEA audit never substitutes nominal P for FDR", {
  df <- data.frame(ID = c("a", "b", "c", "d", "e"), NES = c(1, 2, 3, 4, 5),
                   pvalue = c(0.001, 0.001, 0.001, 0.001, 0.001),
                   p.adjust = c(0.2, 0.05, NA, Inf, -0.01))
  expect_equal(significant_gsea_terms(df)$ID, "b")
  df$p.adjust <- NULL
  expect_equal(nrow(significant_gsea_terms(df)), 0)
  expect_error(audit_gsea_table(df, NA), "finite")
})

test_that("GSEA exports retain unusable rows and explicit zero-significant counts", {
  df <- data.frame(ID = c("known", NA), NES = c(1, NA),
                   pvalue = c(0.001, NA), p.adjust = c(0.1, NA))
  root <- tempfile()
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  file <- file.path(root, "GSEA_GO_test.csv")
  write_gsea_tables(df, file)
  saved <- read.csv(file)
  quality <- read.csv(file.path(root, "GSEA_GO_test_quality.csv"))
  expect_equal(nrow(saved), 2)
  expect_equal(quality$table_rows, 2)
  expect_equal(quality$valid_terms, 1)
  expect_equal(quality$significant_terms, 0)
  expect_identical(saved$significant, c(FALSE, FALSE))
})

test_that("GSEA tests the input-overlap family and BH adjusts that same family", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("gson")
  set.seed(19)
  ranked <- sort(setNames(rnorm(100), as.character(seq_len(100))), decreasing = TRUE)
  sets <- list(large_reference = c(names(ranked)[1:15], paste0("absent", 1:100)),
               valid = names(ranked)[seq(3, 60, by = 6)],
               tiny_overlap = c(names(ranked)[1:2], paste0("absent", 1:30)),
               too_large_overlap = names(ranked)[1:21],
               no_overlap = paste0("absent", 1:30))
  mapping <- do.call(rbind, lapply(names(sets), function(id)
    data.frame(gsid = id, gene = sets[[id]])))
  reference <- gson::gson(mapping,
    gsid2name = data.frame(gsid = names(sets), name = names(sets)),
    species = "synthetic", gsname = "KEGG", keytype = "kegg")
  result <- suppressWarnings(run_kegg_gsea(ranked, organism = reference, min_size = 5,
                                         max_size = 20, seed = 41))
  expect_s4_class(result, "gseaResult")
  df <- as.data.frame(result)
  expect_setequal(df$ID, c("large_reference", "valid"))
  observed <- vapply(sets[df$ID], function(gs) length(intersect(gs, names(ranked))), integer(1))
  expect_equal(df$setSize, unname(observed))
  expect_equal(df$input_overlap, unname(observed))
  expect_equal(df$reference_size[df$ID == "large_reference"], 115)
  expect_equal(df$p.adjust, p.adjust(df$pvalue, "BH"))
  expect_false(any(is.na(df$ID)))
  repeated <- suppressWarnings(run_kegg_gsea(ranked, organism = reference, min_size = 5,
                                           max_size = 20, seed = 41))
  expect_identical(as.data.frame(result), as.data.frame(repeated))
})
