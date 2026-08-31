test_that("TME routing preserves native mouse symbols and maps human references only on request", {
  expr <- matrix(c(10, 20, 30, 40, 15, 25), nrow = 3,
                 dimnames = list(c("Actb", "Gapdh", "Cd274"), c("S1", "S2")))
  native <- prepare_tme_inputs(expr, species = "mouse", need_human = FALSE)
  expect_equal(native$native[rownames(expr), ], as.data.frame(expr))
  expect_null(native$human)
  skip_if_not_installed("babelgene")
  both <- prepare_tme_inputs(expr, species = "mouse", need_human = TRUE)
  expect_equal(both$native, native$native)
  expect_true(all(c("ACTB", "GAPDH", "CD274") %in% rownames(both$human)))
  expect_true(all(c("input_id", "native_symbol", "human_symbol", "source") %in% names(both$mapping)))
})

test_that("TME Ensembl IDs become native symbols before optional ortholog conversion", {
  skip_if_not_installed("org.Mm.eg.db")
  skip_if_not_installed("babelgene")
  ids <- c("ENSMUSG00000029580", "ENSMUSG00000057666")
  expr <- matrix(c(2, 7, 4, 9), 2, dimnames = list(ids, c("A", "B")))
  out <- prepare_tme_inputs(log2(expr + 1), species = "mouse", is_log = TRUE)
  expect_equal(as.numeric(out$native["Actb", ]), as.numeric(expr[1, ]))
  expect_true("ACTB" %in% rownames(out$human))
  expect_equal(out$mapping$input_id, ids)
  annotated <- prepare_tme_inputs(expr, species = "mouse", need_human = FALSE,
                                   native_symbols = setNames(ids, ids))
  expect_equal(annotated$native, out$native)
})

test_that("TME TPM uses the complete feature universe before symbol deduplication", {
  raw <- data.frame(gene_id = c("g1", "g2", "g3"), Length = c(100, 200, 100),
                    S1 = c(10, 20, 1000), S2 = c(20, 40, 1000))
  tpm <- build_tme_tpm(raw, gene_column = "gene_id", count_columns = c("S1", "S2"),
                       length_column = "Length")
  expect_equal(colSums(tpm), c(S1 = 1e6, S2 = 1e6))
  expect_equal(tpm["g1", "S1"], 1e6 * 0.1 / 10.2)
  expect_error(build_tme_tpm(transform(raw, gene_id = c("g1", "g1", "g3")),
                            gene_column = "gene_id", count_columns = c("S1", "S2"),
                            length_column = "Length"), "unique")
})

test_that("TME immune scoring is ssGSEA, with recorded human-signature coverage", {
  skip_if_not_installed("GSVA")
  set.seed(71)
  expr <- matrix(runif(40 * 6, 1, 30), 40,
                 dimnames = list(paste0("G", 1:40), paste0("S", 1:6)))
  sets <- list(A = paste0("g", 1:8), B = paste0("G", 20:30), absent = "MISSING")
  out <- run_tme_ssgsea(expr, gene_sets = sets)
  ref <- GSVA::gsva(GSVA::ssgseaParam(expr, out$gene_sets, minSize = 2,
                                     maxSize = Inf, normalize = TRUE), verbose = FALSE)
  expect_equal(out$scores, ref)
  expect_equal(out$coverage$matched_genes, c(8L, 11L, 0L))
  expect_identical(out$coverage$used, c(TRUE, TRUE, FALSE))
})

test_that("native human and mouse CIBERSORT keep explicit linear scale below 50", {
  for (pkg in c("e1071", "preprocessCore", "future", "furrr", "purrr")) skip_if_not_installed(pkg)
  withr::local_options(list(parallelly.availableCores.methods = "fallback",
                           parallelly.availableCores.fallback = 1L))
  withr::defer(future::plan(future::sequential))
  root <- rprojroot::find_root(rprojroot::is_git_root)
  script <- file.path(root, "references", "CIBERSORT", "CIBERSORT.R")
  for (file in c("LM22.txt", "cibersort_mouse_22.csv")) {
    path <- file.path(root, "references", "CIBERSORT", file)
    sig <- read_cibersort_signature(path)
    weights <- matrix(0, ncol(sig), 2)
    weights[1:2, ] <- matrix(c(.8, .2, .2, .8), 2)
    mix <- sig %*% weights
    colnames(mix) <- c("S1", "S2")
    low <- mix / max(mix) * 5
    high <- low * 100
    a <- run_native_cibersort(low, path, script, perm = 0, QN = FALSE, verbose = FALSE)
    b <- run_native_cibersort(high, path, script, perm = 0, QN = FALSE, verbose = FALSE)
    c <- run_native_cibersort(log2(low + 1), path, script, is_log = TRUE,
                              perm = 0, QN = FALSE, verbose = FALSE)
    expect_equal(a, b, tolerance = 1e-7)
    expect_equal(a, c, tolerance = 1e-7)
    expect_equal(ncol(a) - 4L, ncol(sig))
    expect_equal(unname(rowSums(a[, colnames(sig)])), c(1, 1), tolerance = 1e-7)
    expect_equal(attr(a, "reference_coverage")$matched_genes, nrow(sig))
    expect_true(all(is.na(a[["P-value"]])))
    single <- run_native_cibersort(low[, "S1", drop = FALSE], path, script,
                                    perm = 0, QN = FALSE, verbose = FALSE)
    expect_equal(nrow(single), 1L)
    expect_equal(single$ID, "S1")
    expect_equal(as.numeric(single[1, colnames(sig)]), as.numeric(a[1, colnames(sig)]),
                 tolerance = 1e-7)
    if (file == "cibersort_mouse_22.csv") {
      converted <- low
      rownames(converted) <- toupper(rownames(converted))
      expect_error(run_native_cibersort(converted, path, script, perm = 0, verbose = FALSE),
                   "Insufficient.*overlap")
    }
  }
})

test_that("TME annotation mapping retains raw IDs and duplicate-handling provenance", {
  expr <- matrix(c(4, 20, 8, 6, 30, 9), 3,
                 dimnames = list(c("g1", "g2", "g3"), c("S1", "S2")))
  outdir <- withr::local_tempdir()
  out <- prepare_tme_inputs(expr, species = "mouse", need_human = FALSE, outdir = outdir,
                            native_symbols = c(g1 = "Actb", g2 = "Actb", g3 = "Gapdh"))
  expect_equal(as.numeric(out$native["Actb", ]), as.numeric(expr["g2", ]))
  expect_equal(out$mapping$input_id, rownames(expr))
  expect_equal(out$mapping$source, rep("input_annotation", 3))
  expect_equal(out$mapping$native_retained, c(FALSE, TRUE, TRUE))
  expect_true(all(file.exists(file.path(outdir,
    c("TME_gene_mapping.csv", "TME_input_coverage.csv", "TPM_native_symbols.csv")))))
  expect_false(file.exists(file.path(outdir, "TPM_human_symbols.csv")))
})

test_that("General TME rereads all features and honors the final sample selection", {
  root <- rprojroot::find_root(rprojroot::is_git_root)
  runner <- parse(file = file.path(root, "templates", "General", "run_analysis.R"))
  block <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("if")) &&
                    identical(x[[2]], quote(isTRUE(RUN_TME))) && "tme_outdir" %in% all.names(x),
                    as.list(runner))
  expect_length(block, 1)
  work <- withr::local_tempdir()
  withr::local_dir(work)
  raw <- data.frame(gene_id = c("g1", "g2", "g3"), gene_name = c("ACTB", "GAPDH", "CD274"),
                    gene_biotype = c("protein_coding", "protein_coding", "lncRNA"),
                    Length = c(100, 200, 100), A = c(10, 20, 1000), B = c(20, 30, 2000),
                    C = c(30, 40, 3000))
  write.csv(raw, "input.csv", row.names = FALSE)
  env <- list2env(list(RUN_TME = TRUE, INPUT_FILE = "input.csv", INPUT_FORMAT = "csv",
    GENE_NAME_COL = "gene_name", COUNT_COLS = c("A", "B", "C"),
    SAMPLE_NAMES = c("S1", "S2", "S3"), SPECIES = "human",
    TME_GENE_ID_COLUMN = "gene_id", TME_GENE_LENGTH_COLUMN = "Length", TME_GENE_LENGTH_UNIT = "bp",
    TME_GENE_START_COL = "gene_start", TME_GENE_END_COL = "gene_end", TME_ORTHOLOG_CACHE = NULL,
    RUN_TME_ESTIMATE = FALSE, RUN_TME_IOBR = FALSE, RUN_TME_SSGSEA = FALSE,
    countData = matrix(c(10L, 30L), 1, dimnames = list("ACTB", c("S1", "S3"))),
    colData = data.frame(condition = c("Control", "Treatment"), row.names = c("S1", "S3"))),
    parent = environment())
  eval(block[[1]], env)
  expected <- build_tme_tpm(raw, "gene_id", c("A", "C"), c("S1", "S3"), "Length")
  expect_equal(env$expr_tpm, expected)
  expect_equal(rownames(env$expr_tpm), c("g1", "g2", "g3"))
  expect_equal(colnames(env$expr_tpm), c("S1", "S3"))
  native <- read.csv("4-TME/TPM_native_symbols.csv", row.names = 1, check.names = FALSE)
  expect_equal(native["CD274", "S1"], expected["g3", "S1"])
})

test_that("CIBERSORT comparison requires matching species, signatures and scale settings", {
  sig <- matrix(1:12, 4, dimnames = list(LETTERS[1:4], c("A", "B", "C")))
  expect_true(cibersort_comparison_compatible("human", sig, sig, FALSE, FALSE, input_max = 100))
  expect_false(cibersort_comparison_compatible("mouse", sig, sig, FALSE, FALSE, input_max = 100))
  expect_false(cibersort_comparison_compatible("human", sig, sig * 2, FALSE, FALSE, input_max = 100))
  expect_false(cibersort_comparison_compatible("human", sig, NULL, FALSE, FALSE, input_max = 100))
  expect_false(cibersort_comparison_compatible("human", sig, sig, TRUE, FALSE, input_max = 100))
  expect_false(cibersort_comparison_compatible("human", sig, sig, FALSE, FALSE, input_max = 5))
  expect_false(cibersort_comparison_compatible("human", sig, sig, FALSE, FALSE))
})

test_that("TME runner and notebook feed the native solver species-matched linear input", {
  for (pkg in c("e1071", "preprocessCore", "future", "furrr", "purrr", "jsonlite")) skip_if_not_installed(pkg)
  withr::local_options(list(parallelly.availableCores.methods = "fallback",
                           parallelly.availableCores.fallback = 1L))
  withr::defer(future::plan(future::sequential))
  root <- rprojroot::find_root(rprojroot::is_git_root)
  notebook <- jsonlite::fromJSON(file.path(root, "notebooks", "RNAseq_TME_Deconvolution_Template.ipynb"),
                                  simplifyVector = FALSE)
  notebook_code <- vapply(notebook$cells, function(cell) {
    if (cell$cell_type == "code") paste(unlist(cell$source), collapse = "") else ""
  }, character(1))
  sources <- list(readLines(file.path(root, "templates", "TME", "run_analysis.R")),
                   strsplit(paste(notebook_code, collapse = "\n"), "\n", fixed = TRUE)[[1]])
  work <- withr::local_tempdir()
  meta_path <- file.path(work, "metadata.csv")
  write.csv(data.frame(sample = c("S1", "S2"), condition = c("Control", "Treatment")),
            meta_path, row.names = FALSE)
  for (species in c("human", "mouse")) for (code in sources) {
    sig_path <- file.path(root, "references", "CIBERSORT",
                           if (species == "human") "LM22.txt" else "cibersort_mouse_22.csv")
    sig <- read_cibersort_signature(sig_path)
    mixture <- sig[, 1:2, drop = FALSE]
    mixture <- sweep(mixture, 2, colSums(mixture), "/") * 1e6
    colnames(mixture) <- c("S1", "S2")
    expr_path <- file.path(work, "expression.csv")
    write.csv(data.frame(gene_name = rownames(mixture), log2(mixture + 1), check.names = FALSE),
              expr_path, row.names = FALSE)
    e <- list2env(list(META_FILE = meta_path, SAMPLE_COLUMN = "sample", GROUP_COLUMN = "condition",
      GROUP_LEVELS = NULL, INPUT_MODE = "expression", EXPR_UNIT = "log2_tpm", EXPR_FILE = expr_path,
      GENE_COLUMN = "gene_name", SPECIES = species, ORTHOLOG_CACHE = NULL,
      RUN_ESTIMATE = FALSE, RUN_IOBR = FALSE, RUN_SSGSEA = FALSE,
      OUTDIR = work, TME_DIR = work, CIBERSORT_SIGNATURE = sig_path,
      CIBERSORT_SCRIPT = file.path(root, "references", "CIBERSORT", "CIBERSORT.R"),
      CIBERSORT_PERM = 0, CIBERSORT_QN = FALSE), parent = environment())
    input_block <- parse(text = code[seq.int(grep("^meta <- read_metadata", code)[1],
                                             grep("^expr_tme <- tme_inputs", code)[1])])
    eval(input_block, e)
    expect_equal(as.matrix(e$expr_native[rownames(mixture), , drop = FALSE]), mixture, tolerance = 1e-9)
    expect_null(e$expr_tme)
    start <- grep("^  native_cibersort <- run_native_cibersort", code)[1]
    end <- start - 1L + grep("^  \\)$", code[start:length(code)])[1]
    eval(parse(text = code[start:end]), e)
    expect_equal(e$native_cibersort$sample, c("S1", "S2"))
    expect_equal(attr(e$native_cibersort, "reference_coverage")$matched_genes, nrow(sig))
    expect_equal(unname(rowSums(e$native_cibersort[, colnames(sig)])), c(1, 1), tolerance = 1e-7)
    e$EXPR_UNIT <- "vst"
    expect_error(eval(input_block, e), "VST/rlog is not valid")
    e$INPUT_MODE <- "typo"
    expect_error(eval(input_block, e), "INPUT_MODE must be")
  }
})

test_that("TME cannot declare completion when every requested method produced no result", {
  root <- rprojroot::find_root(rprojroot::is_git_root)
  expressions <- parse(file.path(root, "templates", "TME", "run_analysis.R"))
  block <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("if")) &&
                    "requested_tme_methods" %in% all.names(x), as.list(expressions))
  expect_length(block, 1L)
  e <- list2env(list(requested_tme_methods = c(TRUE, FALSE, FALSE, FALSE),
    estimate_scores = NULL, iobr_results = list(), native_cibersort = NULL,
    ssgsea_scores = NULL), parent = environment())
  expect_error(eval(block[[1]], e), "None of the requested TME methods produced results")
  e$requested_tme_methods[] <- FALSE
  expect_silent(eval(block[[1]], e))
})
