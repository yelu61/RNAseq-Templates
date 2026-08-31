test_that("CLI entry points restore space-mangled Rscript paths", {
  root <- if (exists("repo_root", inherits = TRUE)) {
    get("repo_root", inherits = TRUE)
  } else {
    rprojroot::find_root(rprojroot::is_git_root, path = getwd())
  }
  files <- c(
    file.path("tools", c("build_run_registry.R", "render_report.R")),
    unlist(lapply(c("General", "Limma_Voom", "WGCNA", "TME", "TimeCourse", "TCGA_GEO"),
                  function(topic) file.path("templates", topic,
                                             c("run_analysis.R", "visualize_results.R"))))
  )
  for (file in files) {
    lines <- readLines(file.path(root, file), warn = FALSE)
    expect_true(any(grepl('gsub("~\\\\+~", " ", file_arg)', lines, fixed = TRUE)),
                info = file)
  }
})

test_that("every production CLI rejects old optional outputs before writing anything", {
  root <- if (exists("repo_root", inherits = TRUE)) get("repo_root", inherits = TRUE) else {
    rprojroot::find_root(rprojroot::is_git_root, path = getwd())
  }
  base <- tempfile("protected CLI runs ")
  dir.create(base)
  on.exit(unlink(base, recursive = TRUE), add = TRUE)
  withr::local_envvar(RNASEQ_LIB_DIR = file.path(root, "RNAseq_lib"))
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  for (topic in c("General", "Limma_Voom", "WGCNA", "TME", "TimeCourse", "TCGA_GEO")) {
    run <- file.path(base, topic)
    dir.create(file.path(run, "2-GSEA"), recursive = TRUE)
    writeLines("preserve prior optional result", file.path(run, "2-GSEA", "previous.csv"))
    config <- file.path(run, "config.R")
    # Deliberately omit analysis settings: a stale run must fail before loading
    # dependencies or reading inputs, even when the old module is now disabled.
    writeLines(c("DEFAULT_DEG_PVALUE_COLUMN <- 'padj'", "RUN_TME <- FALSE",
                 paste0("OUTDIR <- ", deparse(if (topic == "WGCNA") file.path(run, "custom_network") else run))),
               config)
    setwd(run)
    before_files <- list.files(run, recursive = TRUE, full.names = TRUE, all.files = TRUE)
    before_md5 <- tools::md5sum(before_files)
    output <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
      c(shQuote(file.path(root, "templates", topic, "run_analysis.R")), shQuote(config)),
      stdout = TRUE, stderr = TRUE))
    expect_identical(attr(output, "status"), 1L, info = topic)
    expect_match(paste(output, collapse = "\n"), "Refusing to reuse", info = topic)
    expect_identical(list.files(run, recursive = TRUE, full.names = TRUE, all.files = TRUE),
                     before_files, info = topic)
    expect_identical(tools::md5sum(before_files), before_md5, info = topic)
  }
})

test_that("production report blocks resolve their template from the backend outside its checkout", {
  root <- if (exists("repo_root", inherits = TRUE)) get("repo_root", inherits = TRUE) else {
    rprojroot::find_root(rprojroot::is_git_root, path = getwd())
  }
  project <- tempfile("external report project ")
  dir.create(project)
  on.exit(unlink(project, recursive = TRUE), add = TRUE)
  withr::local_dir(project)
  expected <- normalizePath(file.path(root, "reports", "analysis_report.qmd"))
  for (topic in c("General", "Limma_Voom", "WGCNA", "TME", "TimeCourse", "TCGA_GEO")) {
    expressions <- parse(file.path(root, "templates", topic, "run_analysis.R"))
    report_blocks <- Filter(function(expr) {
      is.call(expr) && identical(expr[[1]], as.name("if")) &&
        identical(expr[[2]], quote(isTRUE(GENERATE_HTML_REPORT)))
    }, as.list(expressions))
    expect_equal(length(report_blocks), 1L, info = topic)
    env <- new.env(parent = baseenv())
    env$GENERATE_HTML_REPORT <- TRUE
    env$REPORT_TITLE <- "External project"
    env$DEFAULT_THRESHOLD <- "primary"
    env$lib_dir <- file.path(root, "RNAseq_lib")
    env$OUTDIR <- project
    env$RUN_ROOT <- project
    env$requireNamespace <- function(...) TRUE
    env$validate_analysis_report <- function(...) invisible(TRUE)
    captured <- NULL
    env$render_analysis_report <- function(outdir, report_file, template = NULL, params) {
      captured <<- list(template = template, params = params)
      report_file
    }
    # Execute the runner's real report branch and path expressions, replacing
    # only rendering/validation so this regression needs no Quarto or Pandoc.
    invisible(capture.output(eval(report_blocks[[1]], envir = env)))
    expect_identical(captured$template, expected, info = topic)
    expect_true(file.exists(captured$template), info = topic)
    if (topic == "General") expect_identical(captured$params$primary_threshold, "primary")
  }
})
