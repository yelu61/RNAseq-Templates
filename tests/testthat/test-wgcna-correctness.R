# Execute the actual production/notebook analysis blocks on small deterministic
# inputs, rather than maintaining a second implementation in the test.
wgcna_analysis_sources <- function() {
  repo <- rprojroot::find_root(rprojroot::is_git_root)
  paths <- c("templates/WGCNA/run_analysis.R",
             "notebooks/RNAseq_WGCNA_Template.ipynb",
             "examples/demo_RNAseq_WGCNA/RNAseq_WGCNA_Template.ipynb")
  setNames(lapply(paths, function(path) {
    if (grepl("[.]R$", path)) return(readLines(file.path(repo, path), warn = FALSE))
    notebook <- jsonlite::fromJSON(file.path(repo, path), simplifyVector = FALSE)
    code <- vapply(notebook$cells, function(cell) {
      if (cell$cell_type == "code") paste(unlist(cell$source), collapse = "") else ""
    }, character(1))
    strsplit(paste(code, collapse = "\n"), "\n", fixed = TRUE)[[1]]
  }), paths)
}

eval_wgcna_block <- function(code, first, last, envir) {
  begin <- grep(first, code)[1]
  end <- grep(last, code)
  end <- end[end >= begin][1]
  stopifnot(!is.na(begin), !is.na(end))
  eval(parse(text = code[seq.int(begin, end)]), envir = envir)
}

test_that("WGCNA exports one correctly matched eigengene correlation per gene", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("jsonlite")
  set.seed(82)
  datExpr <- matrix(rnorm(16 * 66), nrow = 16,
                    dimnames = list(paste0("S", 1:16), paste0("g", 1:66)))
  # Nonconsecutive labels and a grey gene catch positional/color mismatches.
  labels <- c(rep(2, 60), rep(7, 5), 0)
  net <- list(colors = labels,
              MEs = WGCNA::moduleEigengenes(datExpr, colors = labels)$eigengenes)
  test_dir <- tempfile("wgcna-hubs-")
  dir.create(test_dir)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)
  for (code in wgcna_analysis_sources()) {
    e <- list2env(list(net = net, datExpr = datExpr, TARGET_MODULES = NULL,
                      OUTDIR = tempfile(tmpdir = test_dir)), parent = globalenv())
    for (fun in c("labels2colors", "orderMEs")) e[[fun]] <- getExportedValue("WGCNA", fun)
    for (fun in c("%>%", "arrange", "desc", "bind_rows")) e[[fun]] <- getExportedValue("dplyr", fun)
    dir.create(e$OUTDIR)
    eval_wgcna_block(code, "^moduleColors <-", "^MEs <-", e)
    # Include any explicit eigengene relabelling before the first file export.
    start <- grep("^MEs <-", code)[1] + 1L
    end <- grep("^write.csv\\(data.frame", code)[1] - 1L
    if (start <= end) eval(parse(text = code[start:end]), envir = e)
    expect_setequal(sub("^ME", "", colnames(e$MEs)), unique(e$moduleColors))
    eval_wgcna_block(code, "^gene_module <-", "Hub_genes_all_modules.csv.*row.names", e)
    hubs <- read.csv(file.path(e$OUTDIR, "Hub_genes_all_modules.csv"))
    expect_equal(nrow(hubs), 65L)
    expect_equal(anyDuplicated(hubs$gene), 0L)
    expect_setequal(hubs$gene, colnames(datExpr)[labels != 0])
    for (label in c(2, 7)) {
      color <- WGCNA::labels2colors(label)
      hub <- read.csv(file.path(e$OUTDIR, paste0("Hub_genes_", color, ".csv")))
      expected <- vapply(hub$gene, function(gene) {
        stats::cor(datExpr[, gene], net$MEs[[paste0("ME", label)]])
      }, numeric(1))
      expect_equal(hub$kME, unname(expected), tolerance = 1e-12)
      expect_true(all(diff(abs(hub$kME)) <= 0))
      expect_true(all(hub$module == color))
    }
    for (target in list(c("blue", "grey"), "unmatched")) {
      e$TARGET_MODULES <- target
      e$OUTDIR <- tempfile(tmpdir = test_dir)
      dir.create(e$OUTDIR)
      eval_wgcna_block(code, "^gene_module <-", "Hub_genes_all_modules.csv.*row.names", e)
      selected <- read.csv(file.path(e$OUTDIR, "Hub_genes_all_modules.csv"))
      expect_named(selected, c("gene", "module", "kME"))
      expect_equal(nrow(selected), if ("blue" %in% target) 60L else 0L)
      expect_equal(anyDuplicated(selected$gene), 0L)
      expect_true(all(selected$module == "blue"))
    }
    e$moduleColors <- rep("grey", ncol(datExpr))
    e$TARGET_MODULES <- NULL
    e$OUTDIR <- tempfile(tmpdir = test_dir)
    dir.create(e$OUTDIR)
    eval_wgcna_block(code, "^gene_module <-", "Hub_genes_all_modules.csv.*row.names", e)
    grey <- read.csv(file.path(e$OUTDIR, "Hub_genes_all_modules.csv"))
    expect_named(grey, c("gene", "module", "kME"))
    expect_equal(nrow(grey), 0L)
  }
})

test_that("WGCNA trait p-values use the matching nonmissing sample count", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("jsonlite")
  MEs <- data.frame(MEblue = c(2, 1, 4, 3, 5, 8, 6, 7),
                    MEblack = c(NA, 3, 1, 4, 8, 7, 6, 2))
  traits <- cbind(complete = c(1, 2, 4, 3, 8, 5, 7, 6),
                  partial = c(2, NA, 1, 4, NA, 8, NA, 3),
                  two_only = c(2, 4, rep(NA, 6)))
  test_dir <- tempfile("wgcna-traits-")
  dir.create(test_dir)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)
  for (code in wgcna_analysis_sources()) {
    e <- list2env(list(MEs = MEs, trait_numeric = traits,
                      datExpr = matrix(0, nrow(MEs), 2),
                      corPvalueStudent = WGCNA::corPvalueStudent,
                      OUTDIR = tempfile(tmpdir = test_dir)), parent = globalenv())
    dir.create(e$OUTDIR)
    suppressWarnings(eval_wgcna_block(code, "^moduleTraitCor <-",
      "^write.csv\\(moduleTraitP,", e))
    for (i in seq_len(ncol(MEs))) for (j in 1:2) {
      expected <- stats::cor.test(MEs[[i]], traits[, j])
      expect_equal(unname(e$moduleTraitCor[i, j]), unname(expected$estimate))
      expect_equal(unname(e$moduleTraitP[i, j]), expected$p.value, tolerance = 1e-12)
    }
    expect_true(all(is.na(e$moduleTraitP[, "two_only"])))
    expect_equal(unname(e$moduleTraitN), rbind(c(8, 5, 2), c(7, 4, 1)))
  }
})

test_that("WGCNA re-visualization repairs legacy numeric eigengene hub exports", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dplyr")
  repo <- rprojroot::find_root(rprojroot::is_git_root)
  expressions <- parse(file.path(repo, "templates/WGCNA/visualize_results.R"))
  block <- expressions[vapply(expressions, function(expr) {
    grepl("^if \\(isTRUE\\(DO_HUB\\)\\)", paste(deparse(expr), collapse = " "))
  }, logical(1))]
  stopifnot(length(block) == 1)
  datExpr <- cbind(g1 = c(1, 3, 2, 5, 4), g2 = c(5, 1, 2, 3, 4))
  e <- list2env(list(DO_HUB = TRUE, datExpr = datExpr,
                    MEs = data.frame(ME2 = c(2, 1, 3, 5, 4)),
                    moduleColors = c("blue", "blue"),
                    HUB_MODULES = NULL, HUB_TOP_N = NULL,
                    OUTDIR = tempfile("wgcna-replot-")), parent = globalenv())
  on.exit(unlink(e$OUTDIR, recursive = TRUE), add = TRUE)
  for (fun in c("%>%", "arrange", "desc", "bind_rows")) e[[fun]] <- getExportedValue("dplyr", fun)
  e$labels2colors <- WGCNA::labels2colors
  capture.output(eval(block, envir = e))
  expect_length(e$hub_list, 1)
  expect_equal(e$hub_list$blue$gene, c("g1", "g2"))
  expect_equal(e$hub_list$blue$kME, c(0.7, 0.3))
})
