# Unit tests for data_utils.R

# -----------------------------------------------------------------------------
# read_expression_matrix
# -----------------------------------------------------------------------------

test_that("read_expression_matrix loads CSV and uses first non-numeric column as row names", {
  tmp <- tempfile(fileext = ".csv")
  df <- data.frame(
    gene = c("A", "B", "C"),
    S1 = c(1.0, 2.0, 3.0),
    S2 = c(4.0, 5.0, 6.0),
    stringsAsFactors = FALSE
  )
  utils::write.csv(df, tmp, row.names = FALSE)

  mat <- read_expression_matrix(tmp, file_format = "csv")
  expect_equal(rownames(mat), c("A", "B", "C"))
  expect_equal(colnames(mat), c("S1", "S2"))
  expect_equal(as.numeric(mat["A", ]), c(1, 4))
  unlink(tmp)
})

test_that("read_expression_matrix subsets samples and warns on missing", {
  tmp <- tempfile(fileext = ".csv")
  df <- data.frame(
    gene = c("A", "B"),
    S1 = c(1, 2),
    S2 = c(3, 4),
    stringsAsFactors = FALSE
  )
  utils::write.csv(df, tmp, row.names = FALSE)

  expect_warning(
    mat <- read_expression_matrix(tmp, file_format = "csv", sample_ids = c("S1", "S3")),
    "S3"
  )
  expect_equal(colnames(mat), "S1")
  unlink(tmp)
})

test_that("read_expression_matrix errors on non-numeric expression values", {
  tmp <- tempfile(fileext = ".csv")
  df <- data.frame(
    gene = c("A", "B"),
    S1 = c("1", "x"),
    stringsAsFactors = FALSE
  )
  utils::write.csv(df, tmp, row.names = FALSE)

  expect_error(read_expression_matrix(tmp, file_format = "csv"), "non-numeric")
  unlink(tmp)
})

test_that("read_expression_matrix requires an explicit numeric gene column", {
  tmp <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(entrez_id = c(1, 2), S1 = c(10, 20)), tmp, row.names = FALSE)
  expect_error(read_expression_matrix(tmp, file_format = "csv"), "cannot be inferred safely")
  mat <- read_expression_matrix(tmp, gene_column = "entrez_id", file_format = "csv")
  expect_equal(rownames(mat), c("1", "2"))
  unlink(tmp)
})

test_that("read_expression_matrix deduplicates row names by mean expression", {
  tmp <- tempfile(fileext = ".csv")
  df <- data.frame(
    gene = c("A", "A", "B"),
    S1 = c(1, 10, 5),
    S2 = c(2, 20, 6),
    stringsAsFactors = FALSE
  )
  utils::write.csv(df, tmp, row.names = FALSE)

  mat <- read_expression_matrix(tmp, file_format = "csv")
  expect_equal(rownames(mat), c("A", "B"))
  expect_equal(as.numeric(mat["A", ]), c(10, 20))
  unlink(tmp)
})

# -----------------------------------------------------------------------------
# read_metadata
# -----------------------------------------------------------------------------

test_that("read_metadata validates required columns and factorizes groups", {
  tmp <- tempfile(fileext = ".csv")
  df <- data.frame(
    sample = c("S1", "S2", "S3"),
    condition = c("A", "A", "B"),
    time = c("T0", "T1", "T1"),
    stringsAsFactors = FALSE
  )
  utils::write.csv(df, tmp, row.names = FALSE)

  meta <- read_metadata(
    tmp,
    required_columns = c("condition", "time"),
    group_column = "condition",
    group_levels = c("A", "B"),
    time_column = "time",
    time_levels = c("T0", "T1")
  )
  expect_equal(levels(meta$condition), c("A", "B"))
  expect_equal(levels(meta$time), c("T0", "T1"))
  unlink(tmp)
})

test_that("read_metadata renames duplicated sample columns", {
  tmp <- tempfile(fileext = ".csv")
  df <- data.frame(
    sample = c("S1", "S2"),
    sample = c("S1", "S2"),
    condition = c("A", "B"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  utils::write.csv(df, tmp, row.names = FALSE)

  meta <- read_metadata(tmp)
  expect_true("sample" %in% colnames(meta))
  expect_true("sample_1" %in% colnames(meta))
  expect_equal(meta$sample, c("S1", "S2"))
  unlink(tmp)
})

test_that("read_metadata errors on missing sample column and duplicate IDs", {
  tmp1 <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(id = c("S1", "S2"), condition = c("A", "B")), tmp1, row.names = FALSE)
  expect_error(read_metadata(tmp1), "Sample column")
  unlink(tmp1)

  tmp2 <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(sample = c("S1", "S1"), condition = c("A", "B")), tmp2, row.names = FALSE)
  expect_error(read_metadata(tmp2), "duplicate")
  unlink(tmp2)
})

test_that("read_metadata rejects unrecognized group and time values", {
  tmp <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(sample = c("S1", "S2"), condition = c("A", "C"), time = c("T0", "T2")),
                   tmp, row.names = FALSE)
  expect_error(read_metadata(tmp, group_column = "condition", group_levels = c("A", "B")),
               "not listed in group_levels")
  expect_error(read_metadata(tmp, time_column = "time", time_levels = c("T0", "T1")),
               "not listed in time_levels")
  unlink(tmp)
})

test_that("encode_wgcna_traits excludes identifier mirrors and labels contrasts", {
  traits <- data.frame(
    sample = paste0("S", 1:4),
    sample_1 = paste0("S", 1:4),
    condition = factor(c("Control", "Control", "Treatment", "Treatment"),
                       levels = c("Control", "Treatment")),
    age = c(40, 50, 60, 70),
    invariant = 1,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(traits) <- traits$sample
  out <- encode_wgcna_traits(traits)
  expect_equal(colnames(out), c("condition: Treatment", "age"))
  expect_equal(as.numeric(out[, "condition: Treatment"]), c(0, 0, 1, 1))
  expect_equal(rownames(out), traits$sample)
  expect_type(out, "double")
})

test_that("encode_wgcna_traits returns an empty matrix without varying traits", {
  traits <- data.frame(sample = c("S1", "S2"), batch = c("A", "A"))
  out <- encode_wgcna_traits(traits)
  expect_equal(dim(out), c(2L, 0L))
})

# -----------------------------------------------------------------------------
# validate_samples_match
# -----------------------------------------------------------------------------

test_that("validate_samples_match errors on no overlap and warns on extras", {
  expect_error(
    validate_samples_match(c("S1", "S2"), c("S3", "S4")),
    "No common samples"
  )
  expect_warning(
    validate_samples_match(c("S1", "S2", "S3"), c("S1", "S2", "S4")),
    "S3"
  )
})

test_that("validate_samples_match enforces strict order", {
  expect_error(
    validate_samples_match(c("S1", "S2"), c("S2", "S1"), strict_order = TRUE),
    "order mismatch"
  )
})

# -----------------------------------------------------------------------------
# detect_expression_scale
# -----------------------------------------------------------------------------

test_that("detect_expression_scale classifies raw counts", {
  set.seed(1)
  mat <- matrix(rpois(1000, lambda = 500), nrow = 100)
  res <- detect_expression_scale(mat)
  expect_equal(res$scale, "raw_counts")
  expect_true(res$confidence %in% c("high", "medium"))
})

test_that("detect_expression_scale classifies log2(TPM+1)", {
  mat <- matrix(runif(1000, 0, 12), nrow = 100)
  res <- detect_expression_scale(mat)
  expect_equal(res$scale, "log2_tpm")
})

test_that("detect_expression_scale classifies TPM", {
  # TPM-like: columns sum to ~1e6, no integers, values spread
  mat <- matrix(rexp(1000, rate = 1/100), nrow = 100)
  mat <- t(t(mat) / colSums(mat)) * 1e6
  res <- detect_expression_scale(mat)
  expect_equal(res$scale, "tpm")
})

test_that("detect_expression_scale classifies VST by negative values", {
  mat <- matrix(rnorm(1000, mean = 0, sd = 2), nrow = 100)
  res <- detect_expression_scale(mat)
  expect_equal(res$scale, "vst")
})

# -----------------------------------------------------------------------------
# counts_to_tpm / counts_to_fpkm / extract_gene_lengths
# -----------------------------------------------------------------------------

test_that("counts_to_tpm computes correct TPM values", {
  counts <- matrix(c(100, 200, 300, 400), nrow = 2,
                   dimnames = list(c("A", "B"), c("S1", "S2")))
  lengths <- c(A = 1, B = 2)
  tpm <- counts_to_tpm(counts, lengths)
  expect_equal(colSums(tpm), c(S1 = 1e6, S2 = 1e6))
  expect_true(all(tpm >= 0))
})

test_that("counts_to_tpm rejects zero-expression samples", {
  counts <- matrix(c(0, 0, 10, 20), nrow = 2,
                   dimnames = list(c("A", "B"), c("empty", "S2")))
  expect_error(counts_to_tpm(counts, c(A = 1, B = 2)), "zero total RPK")
})

test_that("counts_to_fpkm computes correct FPKM values", {
  counts <- matrix(c(100, 200, 300, 400), nrow = 2,
                   dimnames = list(c("A", "B"), c("S1", "S2")))
  lengths <- c(A = 1, B = 2)
  fpkm <- counts_to_fpkm(counts, lengths)
  expect_equal(as.numeric(fpkm["A", "S1"]), 100 / 1 / (300 / 1e6))
})

test_that("extract_gene_lengths computes lengths from start/end", {
  annot <- data.frame(
    gene_id = c("A", "B"),
    gene_start = c(1, 1001),
    gene_end = c(1000, 3000),
    stringsAsFactors = FALSE
  )
  lengths <- extract_gene_lengths(annot, id_col = "gene_id")
  expect_equal(lengths, c(A = 1, B = 2))
})

test_that("extract_gene_lengths uses length column in bp", {
  annot <- data.frame(
    gene_id = c("A", "B"),
    length = c(1000, 2000),
    stringsAsFactors = FALSE
  )
  lengths <- extract_gene_lengths(annot, id_col = "gene_id", length_col = "length", length_unit = "bp")
  expect_equal(lengths, c(A = 1, B = 2))
})

# -----------------------------------------------------------------------------
# validate_count_matrix
# -----------------------------------------------------------------------------

test_that("validate_count_matrix passes valid integer counts", {
  mat <- matrix(0:3, nrow = 2, dimnames = list(c("A", "B"), c("S1", "S2")))
  expect_true(validate_count_matrix(mat))
})

test_that("validate_count_matrix catches negatives and non-integers", {
  mat_neg <- matrix(c(0, -1, 2, 3), nrow = 2)
  expect_error(validate_count_matrix(mat_neg), "negative")

  mat_nonint <- matrix(c(0, 1.5, 2, 3), nrow = 2)
  expect_error(validate_count_matrix(mat_nonint), "non-integer")
})

test_that("validate_count_matrix rejects duplicate gene identifiers", {
  mat <- matrix(1:4, nrow = 2, dimnames = list(c("A", "A"), c("S1", "S2")))
  expect_error(validate_count_matrix(mat), "must be unique")
})

# -----------------------------------------------------------------------------
# convert_gene_ids / convert_expression_rownames
# -----------------------------------------------------------------------------

test_that("detect_gene_id_type distinguishes Ensembl and symbols", {
  expect_equal(detect_gene_id_type(c("ENSG000001", "ENSG000002")), "ensembl")
  expect_equal(detect_gene_id_type(c("TP53", "GAPDH")), "symbol")
})

test_that("convert_expression_rownames stops on invalid identifiers", {
  expr <- matrix(1:6, nrow = 3,
                 dimnames = list(c("xxx_not_a_gene_12345", "yyy_not_a_gene_67890", "zzz_not_a_gene_11111"), c("S1", "S2")))
  expect_error(
    convert_expression_rownames(expr, species = "human", target = "symbol"),
    "No gene IDs could be converted"
  )
})

test_that("convert_expression_rownames converts and deduplicates successfully", {
  expr <- matrix(c(1, 10, 2, 2, 20, 4), nrow = 3,
                 dimnames = list(c("id1", "id2", "id3"), c("S1", "S2")))
  original_convert <- get("convert_gene_ids", envir = .GlobalEnv)
  on.exit(assign("convert_gene_ids", original_convert, envir = .GlobalEnv), add = TRUE)
  assign("convert_gene_ids", function(ids, ...) {
    data.frame(id = ids, converted = c("A", "A", "B"), unmapped = FALSE)
  }, envir = .GlobalEnv)
  out <- convert_expression_rownames(expr, species = "human", target = "symbol")
  expect_equal(rownames(out), c("A", "B"))
  expect_equal(as.numeric(out["A", ]), c(10, 20))
})

test_that("validate_expression_contract enforces method input units", {
  tpm <- matrix(c(4e5, 6e5, 3e5, 7e5), nrow = 2,
                dimnames = list(c("A", "B"), c("S1", "S2")))
  expect_silent(validate_expression_contract(tpm, expected = "tpm"))

  counts <- matrix(c(0, 1000, 0, 2000), nrow = 2,
                   dimnames = list(c("A", "B"), c("S1", "S2")))
  expect_error(validate_expression_contract(counts, expected = "vst"), "raw counts")
})

test_that("collapse_by_symbol keeps the max-count row and carries annotation", {
  mat <- matrix(c(1, 1, 1,   10, 20, 30,   2, 2, 2,   5, 5, 5), nrow = 4, byrow = TRUE,
                dimnames = list(NULL, c("S1", "S2", "S3")))
  sym <- c("A", "B", "A", "C")          # A duplicated; B has highest total
  len <- c(100, 2000, 150, 500)
  cc <- collapse_by_symbol(sym, mat, extra = len)
  expect_equal(cc$gene_name, c("B", "C", "A"))   # ordered by total signal desc
  # for duplicated symbol A, the retained row is the one with higher rowSum (row 3: 6 > 3)
  expect_equal(as.numeric(cc$counts[cc$gene_name == "A", ]), c(2, 2, 2))
  expect_equal(cc$extra[cc$gene_name == "A"], 150)  # length carried from the retained A row
  expect_equal(cc$extra[cc$gene_name == "B"], 2000)
  expect_equal(cc$n_in, 4); expect_equal(cc$n_out, 3)
})

test_that("collapse_by_symbol drops NA/empty symbols and works without extra", {
  mat <- matrix(c(1, 2, 3, 4), nrow = 4, ncol = 1)
  sym <- c("A", NA, "", "B")
  cc <- collapse_by_symbol(sym, mat)
  expect_equal(sort(cc$gene_name), c("A", "B"))
  expect_null(cc$extra)
  expect_equal(cc$n_out, 2)
})
