#!/usr/bin/env Rscript
# =============================================================================
# notebook_to_runner.R — convert a template notebook into a config + CLI runner
# =============================================================================
# Turns one of the notebooks/ *.ipynb templates into the headless
# `config.R` + `run_analysis.R` pair (the same shape as templates/General/).
#
# Conversion convention (no papermill tags are used anywhere in this repo):
#   * Parameter cell = the first R code cell after the markdown header
#     "## 1. Parameter Configuration", or the cell carrying the banner comment
#     "# ====... Parameter Configuration ====".
#   * Environment cell = the code cell holding the library()/source(RNAseq_lib)
#     setup; it is folded into the generated runner's bootstrap, not the body.
#   * All remaining code cells are flattened, in order, into the runner body,
#     each preceded by a section banner taken from its preceding markdown cell.
#
# The output is a *draft*: it reproduces the notebook non-interactively, but
# side effects (dir.create), hard-coded OrgDb, display-only calls and OUTDIR
# routing are flagged in the conversion report for a human/agent to refine.
#
# CLI:
#   Rscript tools/notebook_to_runner.R <notebook.ipynb> <out_dir> [--topic NAME]
#
# The functions are also sourceable (e.g. from tests) without running the CLI.
# =============================================================================

options(stringsAsFactors = FALSE)

# ---- Notebook parsing ---------------------------------------------------------

#' Read a .ipynb into a list of cells: list(type, source, meta).
nb_read <- function(path) {
  if (!file.exists(path)) stop("Notebook not found: ", path)
  nb <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  if (is.null(nb$cells)) stop("No `cells` found in notebook: ", path)
  lapply(nb$cells, function(cell) {
    src <- cell$source
    if (is.list(src)) src <- unlist(src, use.names = FALSE)
    if (length(src) == 0) src <- ""
    list(
      type   = cell$cell_type %||% "",
      source = paste(src, collapse = ""),
      meta   = cell$metadata %||% list()
    )
  })
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

.nb_is_code <- function(cell) identical(cell$type, "code")
.nb_is_md   <- function(cell) identical(cell$type, "markdown")

.nb_lines <- function(src) strsplit(src, "\n", fixed = TRUE)[[1]]

# Strip full-line comments so we don't parse commented-out code.
.nb_code_lines <- function(src) {
  lines <- .nb_lines(src)
  lines[!grepl("^\\s*#", lines)]
}

# ---- Cell location ------------------------------------------------------------

#' Locate the parameter cell index (1-based). Errors if not found.
nb_find_param_cell <- function(cells) {
  # 1) markdown header "## 1. Parameter Configuration" -> next code cell.
  md_hit <- which(vapply(cells, function(c)
    .nb_is_md(c) && grepl("#+\\s*1\\s*\\.\\s*Parameter", c$source, ignore.case = TRUE),
    logical(1)))
  if (length(md_hit) > 0) {
    after <- md_hit[1]
    if (after < length(cells)) {
      for (i in seq.int(after + 1L, length(cells))) {
        if (.nb_is_code(cells[[i]])) return(i)
        if (.nb_is_md(cells[[i]]) && grepl("#+\\s*2\\s*\\.", cells[[i]]$source)) break
      }
    }
  }
  # 2) banner comment inside a code cell.
  banner_hit <- which(vapply(cells, function(c)
    .nb_is_code(c) && grepl("=+\\s*Parameter Configuration\\s*=+", c$source, ignore.case = TRUE),
    logical(1)))
  if (length(banner_hit) > 0) return(banner_hit[1])
  # 3) heuristic: a code cell with several top-level assignments + config keywords.
  heur <- which(vapply(cells, function(c) {
    if (!.nb_is_code(c)) return(FALSE)
    code <- .nb_code_lines(c$source)
    n_assign <- sum(vapply(gregexpr("<-", code, fixed = TRUE),
                           function(m) if (m[1] == -1) 0L else length(m), integer(1)))
    n_assign >= 3 && grepl("OUTDIR|INPUT_FILE|SPECIES|EXPR_FILE", c$source)
  }, logical(1)))
  if (length(heur) > 0) return(heur[1])
  stop("Could not locate a parameter cell in the notebook.")
}

#' Locate the environment/setup cell index (library + source(RNAseq_lib)), or NA.
nb_find_env_cell <- function(cells, param_idx) {
  for (i in seq_along(cells)) {
    if (i == param_idx || !.nb_is_code(cells[[i]])) next
    src <- cells[[i]]$source
    has_lib    <- any(grepl("\\blibrary\\(|\\brequire\\(", .nb_code_lines(src)))
    has_source <- grepl("source\\(", src) && grepl("RNAseq_lib|LIB_DIR", src)
    if (has_lib && has_source) return(i)
  }
  # fall back to the code cell right after a "## N. Environment" header
  md_env <- which(vapply(cells, function(c)
    .nb_is_md(c) && grepl("#+\\s*[0-9]+\\s*\\.\\s*Environment", c$source, ignore.case = TRUE),
    logical(1)))
  if (length(md_env) > 0) {
    after <- md_env[1]
    if (after < length(cells)) {
      for (i in seq.int(after + 1L, length(cells))) {
        if (.nb_is_code(cells[[i]])) return(i)
        if (.nb_is_md(cells[[i]])) break
      }
    }
  }
  NA_integer_
}

# ---- Extraction from the environment cell -------------------------------------

#' Package names from library()/require() calls (ignoring commented lines).
nb_extract_libraries <- function(env_src) {
  code <- .nb_code_lines(env_src)
  hits <- regmatches(code, gregexpr("\\b(?:library|require)\\(\\s*([A-Za-z0-9.]+)", code))
  pkgs <- unique(unlist(lapply(hits, function(m)
    if (length(m)) sub(".*\\(\\s*", "", m) else character(0))))
  pkgs[nzchar(pkgs)]
}

#' RNAseq_lib module filenames from source(file.path(LIB_DIR, "x.R")) or source(".../x.R").
nb_extract_lib_modules <- function(env_src) {
  code <- .nb_code_lines(env_src)
  mods <- character(0)
  # source(file.path(LIB_DIR, "plot_utils.R"))
  m1 <- regmatches(code, gregexpr("source\\(file\\.path\\([^,]+,\\s*\"([^\"]+\\.R)\"", code))
  mods <- c(mods, unlist(lapply(m1, function(m)
    if (length(m)) sub(".*,\\s*\"([^\"]+)\".*", "\\1", m) else character(0))))
  # source("path/to/plot_utils.R") or source(file.path(lib_dir, "plot_utils.R"))
  m2 <- regmatches(code, gregexpr("source\\([^)]*\"([^\"/]*\\.R)\"", code))
  mods <- c(mods, unlist(lapply(m2, function(m)
    if (length(m)) sub(".*\"([^\"/]+\\.R)\".*", "\\1", m) else character(0))))
  unique(mods[grepl("\\.R$", mods)])
}

# ---- Config cleaning ----------------------------------------------------------

#' Remove notebook-only side effects from the parameter cell; return list(
#'   config = cleaned source, removed = character vector of removed statements).
nb_clean_config <- function(param_src) {
  lines <- .nb_lines(param_src)
  removed <- character(0)
  keep <- rep(TRUE, length(lines))
  for (i in seq_along(lines)) {
    ln <- lines[i]
    if (grepl("^\\s*rm\\s*\\(\\s*list", ln)) { keep[i] <- FALSE; removed <- c(removed, trimws(ln)) }
    else if (grepl("^\\s*dir\\.create\\s*\\(", ln)) { keep[i] <- FALSE; removed <- c(removed, trimws(ln)) }
    else if (grepl("Configuration (complete|loaded)", ln) && grepl("\\bcat\\(", ln)) { keep[i] <- FALSE; removed <- c(removed, trimws(ln)) }
  }
  list(config = paste(lines[keep], collapse = "\n"), removed = removed)
}

#' Top-level variable names assigned in a config source (for the snapshot block).
nb_config_vars <- function(config_src) {
  code <- .nb_code_lines(config_src)
  lhs <- regmatches(code, regexpr("^\\s*([A-Za-z][A-Za-z0-9._]*)\\s*<-", code))
  vars <- sub("^\\s*([A-Za-z][A-Za-z0-9._]*)\\s*<-.*", "\\1", lhs)
  unique(vars[nzchar(vars)])
}

# ---- Body flattening ----------------------------------------------------------

#' Flatten all code cells except skip_idx into a sectioned body. Returns list(
#'   body = text, sections = data.frame(index, title)).
nb_flatten_body <- function(cells, skip_idx = integer(0)) {
  out <- character(0)
  sections <- data.frame(index = integer(0), title = character(0), stringsAsFactors = FALSE)
  last_title <- ""
  for (i in seq_along(cells)) {
    cell <- cells[[i]]
    if (.nb_is_md(cell)) {
      # capture a clean section title from the markdown header
      hdr <- regmatches(cell$source, regexpr("#+\\s*[^\n]+", cell$source))
      if (length(hdr) && nzchar(hdr)) last_title <- trimws(sub("^#+\\s*", "", hdr))
      next
    }
    if (!.nb_is_code(cell) || i %in% skip_idx) next
    banner_title <- if (nzchar(last_title)) last_title else sprintf("cell %d", i)
    out <- c(out,
      "",
      "# =============================================================================",
      paste0("# ", banner_title),
      "# =============================================================================",
      cell$source)
    sections <- rbind(sections, data.frame(index = i, title = banner_title, stringsAsFactors = FALSE))
  }
  list(body = paste(out, collapse = "\n"), sections = sections)
}

# ---- Lint: suspicious patterns in the flattened body --------------------------

#' Scan body text for patterns that need manual attention in a headless run.
nb_lint_body <- function(body) {
  notes <- character(0)
  add <- function(pattern, msg, perl = FALSE)
    if (grepl(pattern, body, perl = perl)) notes <<- c(notes, msg)
  add("org\\.Hs\\.eg\\.db|org\\.Mm\\.eg\\.db",
      "Hard-coded OrgDb (org.Hs.eg.db / org.Mm.eg.db) in body; parameterize by SPECIES (pick org_db in the runner).")
  add("\\bView\\s*\\(", "Interactive View() call; remove for headless runs.")
  add("\\binstall\\.packages\\s*\\(|BiocManager::install", "Package installation call in body; move to install_dependencies.R.")
  add("(^|[^\\w])setwd\\s*\\(", "setwd() in body changes the production run root; remove it and use configured paths.", perl = TRUE)
  add("\\bdir\\.create\\s*\\(", "dir.create() in body; prefer centralizing output-dir creation in the runner header.")
  add("\\binteractive\\s*\\(", "interactive() branch detected; verify headless behaviour.")
  add("\\bprint\\s*\\(", "print() calls are display-only; harmless but noisy in logs.")
  add("OUTDIR", "Body writes under OUTDIR; set OUTDIR in config.R to route outputs into the numbered run layout.")
  if (length(notes) == 0) notes <- "No obvious headless-unsafe patterns detected."
  notes
}

# ---- Runner assembly ----------------------------------------------------------

.nb_header <- function(topic, notebook) c(
  "#!/usr/bin/env Rscript",
  "# =============================================================================",
  sprintf("# %s — non-interactive analysis pipeline", topic),
  "# =============================================================================",
  sprintf("# Generated by tools/notebook_to_runner.R from: %s", basename(notebook)),
  sprintf("# Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "#",
  "# DRAFT runner: reproduces the notebook headlessly. Review the conversion",
  "# report (conversion_report.txt) for side effects / hard-coded values to refine.",
  "#",
  "#   Rscript run_analysis.R                  # uses ./config.R",
  "#   Rscript run_analysis.R other_config.R   # uses a different config",
  "# =============================================================================",
  "",
  "options(stringsAsFactors = FALSE)"
)

.nb_bootstrap <- function(topic, libraries, modules, config_vars) {
  lib_block <- if (length(libraries))
    paste0("  library(", libraries, ")", collapse = "\n") else "  # (no library() calls detected)"
  src_block <- if (length(modules))
    paste0('for (f in c(', paste(sprintf('"%s"', modules), collapse = ", "), ')) {\n',
           '  source(file.path(lib_dir, f))\n}') else
    "# (no RNAseq_lib source() calls detected; add required modules here)"

  config_objs <- if (length(config_vars))
    paste(sprintf('  "%s"', config_vars), collapse = ",\n") else "  character(0)"

  c(
    "",
    "# ---- Locate this script and the config file ----------------------------------",
    "invocation_dir <- normalizePath(getwd(), mustWork = TRUE)",
    "cmd_args <- commandArgs(trailingOnly = FALSE)",
    'file_arg <- sub("^--file=", "", grep("^--file=", cmd_args, value = TRUE))',
    "script_dir <- if (length(file_arg) > 0 && nzchar(file_arg)) dirname(normalizePath(file_arg[1])) else invocation_dir",
    "",
    "user_args <- commandArgs(trailingOnly = TRUE)",
    'config_path <- if (length(user_args) >= 1) normalizePath(path.expand(user_args[1]), mustWork = FALSE) else file.path(script_dir, "config.R")',
    "if (!file.exists(config_path)) {",
    '  stop("Config file not found: ", config_path,',
    '       "\\nProvide one as: Rscript run_analysis.R path/to/config.R")',
    "}",
    "config_path <- normalizePath(config_path, mustWork = TRUE)",
    "",
    'cat("========================================\\n")',
    sprintf('cat("%s — run_analysis.R\\n")', topic),
    'cat("Working dir :", getwd(), "\\n")',
    'cat("Config file :", config_path, "\\n")',
    'cat("========================================\\n\\n")',
    "",
    "source(config_path, local = globalenv())",
    "",
    "# ---- Resolve RNAseq_lib -------------------------------------------------------",
    'lib_dir <- Sys.getenv("RNASEQ_LIB_DIR", unset = NA_character_)',
    "if (is.na(lib_dir) || !dir.exists(lib_dir)) {",
    '  if (dir.exists(file.path(invocation_dir, "RNAseq_lib"))) {',
    '    lib_dir <- file.path(invocation_dir, "RNAseq_lib")',
    '  } else if (dir.exists(file.path(script_dir, "RNAseq_lib"))) {',
    '    lib_dir <- file.path(script_dir, "RNAseq_lib")',
    "  } else {",
    "    repo_root <- tryCatch(rprojroot::find_root(rprojroot::is_git_root, path = script_dir), error = function(e) NA_character_)",
    '    if (!is.na(repo_root) && dir.exists(file.path(repo_root, "RNAseq_lib"))) {',
    '      lib_dir <- file.path(repo_root, "RNAseq_lib")',
    '    } else if (dir.exists(file.path(dirname(script_dir), "RNAseq_lib"))) {',
    '      lib_dir <- file.path(dirname(script_dir), "RNAseq_lib")',
    "    } else {",
    '      stop("Could not locate RNAseq_lib. Set RNASEQ_LIB_DIR or run alongside RNAseq_lib/.")',
    "    }",
    "  }",
    "}",
    'cat("RNAseq_lib  :", normalizePath(lib_dir), "\\n\\n")',
    "",
    "# ---- Load libraries -----------------------------------------------------------",
    "suppressPackageStartupMessages({",
    lib_block,
    "})",
    "",
    "# ---- Source RNAseq_lib modules ------------------------------------------------",
    src_block,
    "theme_set(theme_publication())",
    "",
    "# ---- Config snapshot ------------------------------------------------------------",
    'dir.create("0-Config", showWarnings = FALSE, recursive = TRUE)',
    "config_objects <- c(",
    config_objs,
    ")",
    "config_lines <- c(",
    sprintf('  "# %s analysis configuration snapshot",', topic),
    '  paste0("# Saved: ", Sys.time()),',
    "  unlist(lapply(config_objects, function(obj) {",
    '    if (exists(obj, inherits = TRUE)) c(paste0("\\n", obj, " <- "), capture.output(dput(get(obj, inherits = TRUE)))) else NULL',
    "  }))",
    ")",
    'writeLines(config_lines, "./0-Config/analysis_config_used.R")'
  )
}

#' Assemble the full run_analysis.R text.
nb_build_runner <- function(topic, notebook, libraries, modules, config_vars, body) {
  paste(c(
    .nb_header(topic, notebook),
    .nb_bootstrap(topic, libraries, modules, config_vars),
    "",
    "# =============================================================================",
    "# Pipeline body (flattened from notebook code cells)",
    "# =============================================================================",
    body,
    "",
    'cat("\\n========================================\\n")',
    sprintf('cat("%s analysis COMPLETE.\\n")', topic),
    'cat("========================================\\n")',
    ""
  ), collapse = "\n")
}

#' Assemble config.R text with a standard header.
nb_build_config <- function(topic, notebook, cleaned_config) {
  paste(c(
    "# =============================================================================",
    sprintf("# %s — Analysis Configuration", topic),
    "# =============================================================================",
    "# Generated by tools/notebook_to_runner.R from the notebook parameter cell:",
    sprintf("#   %s", basename(notebook)),
    "# This is the ONLY file you need to edit for a new project. Review the",
    "# removed side effects listed in conversion_report.txt (dir.create etc.) —",
    "# output directories are created by run_analysis.R instead.",
    "# =============================================================================",
    "",
    "options(stringsAsFactors = FALSE)",
    "",
    cleaned_config,
    ""
  ), collapse = "\n")
}

# ---- Top-level conversion -----------------------------------------------------

#' Convert one notebook -> out_dir/{config.R, run_analysis.R, conversion_report.txt}.
#' Returns (invisibly) a list describing the conversion.
nb_convert <- function(notebook, out_dir, topic = NULL) {
  cells <- nb_read(notebook)
  if (is.null(topic) || !nzchar(topic)) {
    topic <- tools::file_path_sans_ext(basename(notebook))
    topic <- gsub("^RNAseq_", "", topic)
    topic <- gsub("_Template$", "", topic)
  }

  param_idx <- nb_find_param_cell(cells)
  env_idx   <- nb_find_env_cell(cells, param_idx)
  skip_idx  <- c(param_idx, if (!is.na(env_idx)) env_idx)

  env_src <- if (!is.na(env_idx)) cells[[env_idx]]$source else ""
  libraries <- nb_extract_libraries(env_src)
  modules   <- nb_extract_lib_modules(env_src)

  cleaned  <- nb_clean_config(cells[[param_idx]]$source)
  cfg_vars <- nb_config_vars(cleaned$config)

  flat <- nb_flatten_body(cells, skip_idx = skip_idx)
  lint <- nb_lint_body(flat$body)

  config_text <- nb_build_config(topic, notebook, cleaned$config)
  runner_text <- nb_build_runner(topic, notebook, libraries, modules, cfg_vars, flat$body)

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  writeLines(config_text, file.path(out_dir, "config.R"))
  writeLines(runner_text, file.path(out_dir, "run_analysis.R"))

  report <- c(
    "==========================================================",
    sprintf("Conversion report — %s", topic),
    "==========================================================",
    sprintf("Source notebook : %s", normalizePath(notebook)),
    sprintf("Output dir      : %s", normalizePath(out_dir)),
    sprintf("Generated       : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    sprintf("Parameter cell  : cell %d", param_idx),
    sprintf("Environment cell: %s", if (is.na(env_idx)) "(none detected)" else sprintf("cell %d", env_idx)),
    sprintf("Libraries       : %s", if (length(libraries)) paste(libraries, collapse = ", ") else "(none)"),
    sprintf("Lib modules     : %s", if (length(modules)) paste(modules, collapse = ", ") else "(none)"),
    sprintf("Config variables: %d", length(cfg_vars)),
    "",
    "-- Removed notebook side effects (now handled by the runner) --",
    if (length(cleaned$removed)) paste0("  * ", cleaned$removed) else "  (none)",
    "",
    "-- Flattened body sections --",
    if (nrow(flat$sections)) sprintf("  [%d] %s", flat$sections$index, flat$sections$title) else "  (none)",
    "",
    "-- Patterns to review before production use --",
    paste0("  ! ", lint),
    "",
    "Next steps:",
    "  1. Review config.R (set real paths/samples/groups).",
    "  2. Review run_analysis.R body against the conversion notes above.",
    "  3. Route outputs into the numbered run layout (0-Config/1-DEG/...).",
    "  4. Add Analysis_summary.txt / sessionInfo.txt / optional HTML report tail.",
    ""
  )
  writeLines(report, file.path(out_dir, "conversion_report.txt"))

  invisible(list(
    topic = topic, out_dir = out_dir, param_cell = param_idx, env_cell = env_idx,
    libraries = libraries, modules = modules, config_vars = cfg_vars,
    removed = cleaned$removed, sections = flat$sections, lint = lint
  ))
}

# ---- CLI ----------------------------------------------------------------------

.nb_cli <- function(argv) {
  usage <- paste(
    "Usage: Rscript tools/notebook_to_runner.R <notebook.ipynb> <out_dir> [--topic NAME]",
    sep = "\n")
  if (length(argv) < 2) stop(usage, call. = FALSE)
  notebook <- argv[1]
  out_dir  <- argv[2]
  topic    <- NULL
  ti <- which(argv == "--topic")
  if (length(ti) > 0 && length(argv) >= ti + 1) topic <- argv[ti + 1]

  res <- nb_convert(notebook, out_dir, topic = topic)
  cat("Converted '", basename(notebook), "' -> ", normalizePath(out_dir), "\n", sep = "")
  cat("  topic        :", res$topic, "\n")
  cat("  param cell   :", res$param_cell, "\n")
  cat("  config vars  :", length(res$config_vars), "\n")
  cat("  body sections:", nrow(res$sections), "\n")
  cat("  review notes :", length(res$lint), " (see conversion_report.txt)\n")
  invisible(res)
}

# Run the CLI only when executed as a script (not when source()d by tests).
if (sys.nframe() == 0L) {
  .nb_cli(commandArgs(trailingOnly = TRUE))
}
