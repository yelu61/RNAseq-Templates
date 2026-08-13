# Tests for tools/notebook_to_runner.R — the notebook -> config + runner converter.
#
# The converter is sourceable without running its CLI (the CLI only fires when
# sys.nframe() == 0). These tests build small synthetic notebooks in temp files
# and exercise the parsing / extraction / flattening / lint / end-to-end paths.

repo_root <- tryCatch(
  rprojroot::find_root(rprojroot::is_git_root),
  error = function(e) getwd()
)
source(file.path(repo_root, "tools", "notebook_to_runner.R"))

# --- helpers ------------------------------------------------------------------

# Write a minimal .ipynb with the given cells (list of list(type, source)).
write_nb <- function(cells) {
  path <- tempfile(fileext = ".ipynb")
  json_cells <- lapply(cells, function(c) {
    if (c$type == "code") {
      list(cell_type = "code", source = c$source, metadata = list(),
           outputs = list(), execution_count = NULL)
    } else {
      list(cell_type = "markdown", source = c$source, metadata = list())
    }
  })
  jsonlite::write_json(
    list(cells = json_cells, metadata = list(), nbformat = 4L, nbformat_minor = 5L),
    path, auto_unbox = TRUE, pretty = TRUE)
  path
}

md <- function(txt) list(type = "markdown", source = txt)
code <- function(txt) list(type = "code", source = txt)

# A canonical template-shaped notebook: header, param cell, environment, body.
make_template_nb <- function() {
  write_nb(list(
    md("# Demo Template"),
    md("## 1. Parameter Configuration"),
    code(paste(
      "# ===================== Parameter Configuration =====================",
      "SPECIES <- \"human\"",
      "INPUT_FILE <- \"./0-Data/counts.tsv\"",
      "OUTDIR <- \"Demo_Output\"",
      "dir.create(OUTDIR, showWarnings = FALSE)",
      "cat(\"Configuration complete!\")", sep = "\n")),
    md("## 2. Environment"),
    code(paste(
      "suppressPackageStartupMessages({ library(edgeR); library(limma) })",
      "LIB_DIR <- \"../RNAseq_lib\"",
      "source(file.path(LIB_DIR, \"plot_utils.R\"))",
      "source(file.path(LIB_DIR, \"deg_utils.R\"))", sep = "\n")),
    md("## 3. Load Data"),
    code("rawcount <- read_count_table(INPUT_FILE)"),
    md("## 4. Run"),
    code("print(head(rawcount))")
  ))
}

# --- nb_read ------------------------------------------------------------------

test_that("nb_read parses cells into type/source", {
  cells <- nb_read(make_template_nb())
  expect_equal(length(cells), 9)
  expect_equal(cells[[1]]$type, "markdown")
  expect_equal(cells[[3]]$type, "code")
  expect_match(cells[[3]]$source, "SPECIES")
})

test_that("nb_read errors on a missing file", {
  expect_error(nb_read("/nonexistent/path.ipynb"), "Notebook not found")
})

# --- parameter cell location ---------------------------------------------------

test_that("nb_find_param_cell uses the '## 1. Parameter Configuration' header", {
  cells <- nb_read(make_template_nb())
  expect_equal(nb_find_param_cell(cells), 3)
})

test_that("nb_find_param_cell falls back to the banner comment", {
  nb <- write_nb(list(
    md("# No numbered header here"),
    code("# ==== Parameter Configuration ====\nSPECIES <- \"mouse\"\nX <- 1")
  ))
  expect_equal(nb_find_param_cell(nb_read(nb)), 2)
})

test_that("nb_find_param_cell falls back to the assignment heuristic", {
  nb <- write_nb(list(
    code("library(edgeR)"),  # not a param cell
    code("INPUT_FILE <- \"a\"\nSPECIES <- \"human\"\nOUTDIR <- \"o\"\nGROUPS <- c(\"a\",\"b\")")
  ))
  expect_equal(nb_find_param_cell(nb_read(nb)), 2)
})

test_that("nb_find_param_cell errors when nothing matches", {
  nb <- write_nb(list(code("x <- 1"), code("print(x)")))
  expect_error(nb_find_param_cell(nb_read(nb)), "Could not locate")
})

test_that("nb_find_param_cell handles a terminal parameter header without indexing past the notebook", {
  cells <- list(md("# Demo"), md("## 1. Parameter Configuration"))
  expect_error(nb_find_param_cell(cells), "Could not locate")
})

# --- environment cell + extraction --------------------------------------------

test_that("nb_find_env_cell finds the library+source setup cell", {
  cells <- nb_read(make_template_nb())
  expect_equal(nb_find_env_cell(cells, param_idx = 3), 5)
})

test_that("nb_find_env_cell handles a terminal environment header", {
  cells <- list(code("SPECIES <- 'human'"), md("## 2. Environment"))
  expect_true(is.na(nb_find_env_cell(cells, param_idx = 1)))
})

test_that("nb_extract_libraries ignores commented lines", {
  env <- "# library(commented)\nsuppressPackageStartupMessages({ library(edgeR)\nlibrary(limma) })"
  expect_setequal(nb_extract_libraries(env), c("edgeR", "limma"))
})

test_that("nb_extract_lib_modules pulls module filenames", {
  env <- "source(file.path(LIB_DIR, \"plot_utils.R\"))\nsource(file.path(lib_dir, \"deg_utils.R\"))"
  expect_setequal(nb_extract_lib_modules(env), c("plot_utils.R", "deg_utils.R"))
})

# --- config cleaning -----------------------------------------------------------

test_that("nb_clean_config strips rm/dir.create/completion-cat side effects", {
  src <- paste(
    "rm(list = ls())",
    "SPECIES <- \"human\"",
    "OUTDIR <- \"o\"",
    "dir.create(OUTDIR)",
    "cat(\"Configuration complete!\")", sep = "\n")
  out <- nb_clean_config(src)
  expect_match(out$config, "SPECIES")
  expect_false(grepl("rm\\(list", out$config))
  expect_false(grepl("dir\\.create", out$config))
  expect_false(grepl("Configuration complete", out$config))
  expect_true(any(grepl("dir\\.create", out$removed)))
})

test_that("nb_config_vars extracts top-level assignments", {
  src <- "SPECIES <- \"human\"\nif (TRUE) {\n  nested <- 1\n}\nCOMPARISONS <- list(c(\"a\",\"b\",\"c\"))"
  vars <- nb_config_vars(src)
  expect_true(all(c("SPECIES", "COMPARISONS") %in% vars))
})

# --- body flattening -----------------------------------------------------------

test_that("nb_flatten_body skips param+env and inserts section banners", {
  cells <- nb_read(make_template_nb())
  flat <- nb_flatten_body(cells, skip_idx = c(3, 5))
  expect_match(flat$body, "read_count_table")          # body cell kept
  expect_false(grepl("SPECIES <-", flat$body))          # param cell skipped
  expect_false(grepl("suppressPackageStartupMessages", flat$body))  # env cell skipped
  expect_match(flat$body, "# 3. Load Data")             # section banner from markdown
  expect_equal(nrow(flat$sections), 2)                  # two body code cells
})

# --- lint ----------------------------------------------------------------------

test_that("nb_lint_body flags hard-coded OrgDb and OUTDIR", {
  notes <- nb_lint_body("res <- map_symbols_to_entrez(x, org.Hs.eg.db)\nwrite.csv(y, file.path(OUTDIR, \"x.csv\"))")
  expect_true(any(grepl("OrgDb", notes)))
  expect_true(any(grepl("OUTDIR", notes)))
})

test_that("nb_lint_body returns an all-clear note for clean bodies", {
  notes <- nb_lint_body("x <- 1\ny <- x + 1")
  expect_equal(notes, "No obvious headless-unsafe patterns detected.")
})

# --- end-to-end conversion -----------------------------------------------------

test_that("nb_convert writes config.R, run_analysis.R and a conversion report", {
  out <- file.path(tempdir(), paste0("nbconv_", as.integer(runif(1, 1, 1e9))))
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  res <- nb_convert(make_template_nb(), out, topic = "Demo")

  expect_true(file.exists(file.path(out, "config.R")))
  expect_true(file.exists(file.path(out, "run_analysis.R")))
  expect_true(file.exists(file.path(out, "conversion_report.txt")))

  cfg <- readLines(file.path(out, "config.R"), warn = FALSE)
  expect_true(any(grepl("SPECIES", cfg)))
  # no executable dir.create (header prose may mention it in a comment)
  cfg_code <- cfg[!grepl("^\\s*#", cfg)]
  expect_false(any(grepl("dir\\.create", cfg_code)))  # side effect moved out

  runner <- readLines(file.path(out, "run_analysis.R"), warn = FALSE)
  expect_true(any(grepl("RNASEQ_LIB_DIR", runner)))
  expect_true(any(grepl("invocation_dir", runner)))
  expect_false(any(grepl("setwd\\(script_dir\\)", runner)))
  expect_true(any(grepl("analysis_config_used\\.R", runner)))
  expect_true(any(grepl("library\\(edgeR\\)", runner)))     # libraries carried into bootstrap
  expect_true(any(grepl("read_count_table", runner)))       # body carried over
  expect_equal(res$param_cell, 3)
  expect_equal(res$env_cell, 5)
})

test_that("nb_convert runs on the shipped limma-voom notebook", {
  nb <- file.path(repo_root, "notebooks", "RNAseq_limma_voom_Template.ipynb")
  skip_if_not(file.exists(nb), "limma notebook not present")
  out <- file.path(tempdir(), paste0("nbconv_limma_", as.integer(runif(1, 1, 1e9))))
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  res <- nb_convert(nb, out)
  expect_equal(res$param_cell, 3)
  expect_true("limma_voom_utils.R" %in% res$modules)
  expect_true(file.exists(file.path(out, "config.R")))
})

test_that("all production runners preserve the invocation directory and resolve configs first", {
  runner_dirs <- c("General", "Limma_Voom", "TCGA_GEO", "TME", "TimeCourse", "WGCNA")
  for (topic in runner_dirs) {
    path <- file.path(repo_root, "templates", topic, "run_analysis.R")
    expect_true(file.exists(path), info = topic)
    runner <- readLines(path, warn = FALSE)
    expect_true(any(grepl("invocation_dir <- normalizePath\\(getwd", runner)), info = topic)
    expect_true(any(grepl("normalizePath\\(path.expand\\(user_args\\[1\\]", runner)), info = topic)
    expect_false(any(grepl("setwd\\(script_dir\\)", runner)), info = topic)
  }
})
