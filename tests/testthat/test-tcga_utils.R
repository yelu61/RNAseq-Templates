# Tests for tcga_utils.R (pure helpers; GDC download not tested)

test_that("infer_tcga_tumor_normal classifies by sample-type code", {
  barcodes <- c(
    "TCGA-XX-0001-01A",  # 01 -> tumor
    "TCGA-XX-0002-11A",  # 11 -> normal
    "TCGA-XX-0003-02B",  # 02 -> tumor
    "TCGA-XX-0004-10A"   # 10 -> normal
  )
  out <- infer_tcga_tumor_normal(barcodes)
  expect_equal(out, c("Tumor", "Normal", "Tumor", "Normal"))

  custom <- infer_tcga_tumor_normal(barcodes, tumor_label = "T", normal_label = "N")
  expect_equal(custom, c("T", "N", "T", "N"))
})

test_that("symbolize_and_dedup keeps highest-expression symbol per group", {
  skip_if_not_installed("dplyr")
  # rows: ENSG1/ENSG2 -> GeneA ; ENSG3/ENSG4 -> GeneB (S1,S2 columns)
  mat <- matrix(c(1, 5, 2, 9, 3, 7, 4, 8), nrow = 4, byrow = TRUE,
                dimnames = list(c("ENSG1", "ENSG2", "ENSG3", "ENSG4"), c("S1", "S2")))
  id_map <- data.frame(
    gene_id = c("ENSG1", "ENSG2", "ENSG3", "ENSG4"),
    gene_name = c("GeneA", "GeneA", "GeneB", "GeneB"),
    stringsAsFactors = FALSE
  )
  out <- symbolize_and_dedup(mat, id_map)
  expect_equal(sort(rownames(out)), c("GeneA", "GeneB"))
  expect_equal(colnames(out), c("S1", "S2"))
  # GeneA keeps ENSG2 (mean (2+9)/2=5.5) over ENSG1 ((1+5)/2=3): S1 == 2
  expect_equal(unname(out["GeneA", "S1"]), 2)
  expect_equal(unname(out["GeneA", "S2"]), 9)
  # GeneB keeps ENSG4 (mean (4+8)/2=6) over ENSG3 ((3+7)/2=5): S1 == 4
  expect_equal(unname(out["GeneB", "S1"]), 4)
  expect_equal(unname(out["GeneB", "S2"]), 8)
})

test_that("build_id_map_from_se extracts and cleans the id->symbol map", {
  skip_if_not_installed("SummarizedExperiment")
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:6, nrow = 3, dimnames = list(c("E1", "E2", "E3"), c("S1", "S2")))),
    rowData = S4Vectors::DataFrame(
      gene_id = c("E1", "E2", "E3"),
      gene_name = c("GeneA", "", "GeneC"),
      gene_type = c("protein_coding", "protein_coding", "lncRNA")
    )
  )
  id_map <- build_id_map_from_se(se)
  # E2 dropped because its symbol is empty.
  expect_equal(nrow(id_map), 2)
  expect_true(all(c("gene_id", "gene_name", "gene_type") %in% colnames(id_map)))
  expect_false("" %in% id_map$gene_name)
})

test_that("infer_tcga_tumor_normal handles code boundary at 10", {
  # Use full-length barcodes so substr(., 14, 15) captures the sample-type code.
  expect_equal(infer_tcga_tumor_normal("TCGA-AB-1234-09A"), "Tumor")   # 09 < 10
  expect_equal(infer_tcga_tumor_normal("TCGA-AB-1234-10A"), "Normal")  # 10 not < 10
  # Barcodes too short to carry a code yield NA (no spurious classification).
  expect_true(is.na(infer_tcga_tumor_normal("TCGA-A-B-09X")))
})
