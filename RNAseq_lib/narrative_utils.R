# Narrative-report helpers for RNAseq-Templates.
#
# Backfilled from the bespoke Chinese narrative reports of 202607XXR_RPE and
# 202605CF (reports/CF_analysis_report.Rmd). A narrative report is a
# project-specific, conclusion-led document (see references/NARRATIVE_REPORT.md)
# layered on top of the generic technical report. These helpers provide the
# shared machinery so each project only writes prose, not plumbing:
#
#   1. make_report_citers() — inline-citation closures (sv/gl/gp/gnes/ora_terms)
#      that read numbers from saved run outputs at render time. The report
#      contract forbids hand-typed headline values; always cite through these.
#   2. show_pdf() — embed a run-bundle PDF figure (first page) as base64 PNG.
#   3. validate_narrative_report() — the QA gate a rendered narrative report
#      must pass before it is curated into results/reports/.
#   4. scaffold_narrative_report() — copy the generic scaffold Rmd into a
#      project so a new narrative report starts from the approved skeleton.

# ---------------------------------------------------------------------------
# Inline citers
# ---------------------------------------------------------------------------

# Build the citation closures for one run bundle. `run_dir` must follow the
# standard template layout (1-DEG/, 2-GSEA/, 3-Visualization/). Returns a named
# list; in a report chunk do `list2env(make_report_citers(run_dir), environment())`
# to bring sv/gl/gp/gnes/ora_terms/show_pdf/show_table into scope for all chunks.
make_report_citers <- function(run_dir) {
  read_one <- function(...) {
    path <- file.path(run_dir, ...)
    if (!file.exists(path)) return(NULL)
    tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
             error = function(e) NULL)
  }
  deg_summary <- read_one("1-DEG", "DEG_threshold_summary.csv")
  all_genes_cache <- new.env(parent = emptyenv())

  load_all_genes <- function(comp) {
    if (exists(comp, envir = all_genes_cache)) return(get(comp, envir = all_genes_cache))
    x <- read_one("1-DEG", "all_genes", paste0("DESeq2_all_genes_", comp, ".csv"))
    assign(comp, x, envir = all_genes_cache)
    x
  }
  fmt_p <- function(x) ifelse(x < 0.001, formatC(x, format = "e", digits = 2),
                              formatC(x, format = "f", digits = 3))
  # Gene-level value (default shrunken log2FC) from the all-gene table.
  deg_val <- function(gene, comp, col = "log2FoldChange_shrunken") {
    x <- load_all_genes(comp)
    if (is.null(x) || !col %in% colnames(x)) return(NA_real_)
    v <- x[[col]][x$gene_name == gene]
    if (!length(v) || is.na(v[1])) NA_real_ else as.numeric(v[1])
  }

  list(
    # DEG count for one threshold/contrast cell of the threshold grid.
    sv = function(threshold, comp, col = "Total") {
      if (is.null(deg_summary)) return("NA")
      v <- deg_summary[deg_summary$Threshold == threshold & deg_summary$Comparison == comp, col]
      if (!length(v)) return("NA")
      as.character(v[[1]])
    },
    deg_val = deg_val,
    gl = function(gene, comp) sprintf("%+.1f", deg_val(gene, comp)),
    gp = function(gene, comp) {
      v <- deg_val(gene, comp, "padj")
      if (is.na(v)) "NA" else fmt_p(v)
    },
    # GSEA NES of the first term whose Description contains `pattern` (fixed).
    gnes = function(pattern, comp, db = c("GO", "KEGG")) {
      db <- match.arg(db)
      x <- read_one("2-GSEA", paste0("GSEA_", db, "_", comp, ".csv"))
      if (is.null(x) || !all(c("Description", "NES") %in% colnames(x))) return("NA")
      hit <- x$NES[grepl(pattern, x$Description, fixed = TRUE)]
      if (!length(hit) || is.na(hit[1])) "NA" else sprintf("%+.2f", as.numeric(hit[1]))
    },
    # Top ORA term labels for one direction/contrast at a given threshold tier.
    ora_terms = function(comp, direction = c("UP", "DOWN"), db = c("GO", "KEGG"),
                         n = 5, threshold = "standard") {
      direction <- match.arg(direction); db <- match.arg(db)
      x <- read_one("2-GSEA", threshold, paste0(db, "_ORA_", direction, "_", comp, ".csv"))
      if (is.null(x) || !"Description" %in% colnames(x) || !nrow(x)) return("（无显著条目）")
      x <- x[order(x$p.adjust), ]
      paste(utils::head(x$Description, n), collapse = "；")
    }
  )
}

# ---------------------------------------------------------------------------
# Figure embedding
# ---------------------------------------------------------------------------

# Render the first PDF page to PNG and embed it as a base64 <figure>. Browser
# PDF objects are unreliable; degrades to a text note when the figure or
# pdftoppm is missing. `base_dir` is prepended to `rel` (usually the run dir).
show_pdf <- function(rel, caption, base_dir = ".", width = "100%",
                     missing_note = "_图不可用：") {
  full <- file.path(base_dir, rel)
  if (!file.exists(full)) {
    cat(missing_note, rel, "_\n", sep = "")
    return(invisible(NULL))
  }
  if (!nzchar(Sys.which("pdftoppm"))) {
    cat("_PDF 预览不可用（未安装 pdftoppm）：", rel, "_\n", sep = "")
    return(invisible(NULL))
  }
  base <- tempfile(pattern = "narrative-report-")
  png <- paste0(base, ".png")
  code <- suppressWarnings(system2(
    "pdftoppm",
    c("-f", "1", "-l", "1", "-singlefile", "-r", "150", "-png",
      shQuote(full), shQuote(base)),
    stdout = FALSE, stderr = FALSE
  ))
  if (!identical(code, 0L) || !file.exists(png)) {
    cat("_PDF 预览失败：", rel, "_\n", sep = "")
    return(invisible(NULL))
  }
  uri <- knitr::image_uri(png)
  cat(sprintf(
    '<figure style="margin:1.2em 0"><img src="%s" style="width:%s;height:auto;border:1px solid #e5e5e5"><figcaption style="color:#555;font-size:0.92em;margin-top:0.5em">%s</figcaption></figure>',
    uri, width, caption
  ))
}

# ---------------------------------------------------------------------------
# QA gate
# ---------------------------------------------------------------------------

# Validate a rendered narrative HTML report before curation. Hard failures
# (stop): file missing, fewer than `min_figures` embedded figures, explicit
# missing-figure markers, missing required section headings. Suspicious inline
# "NA" values (a failed citer lookup, e.g. "NES ≈ NA") are reported as
# warnings — they can be legitimate (a gene absent from a table) but must be
# eyeballed. Returns invisibly a list with counts, invisibly TRUE on pass.
validate_narrative_report <- function(html_path,
                                      min_figures = 1,
                                      required_sections = character(0),
                                      missing_markers = c("图不可用", "PDF 预览失败",
                                                          "PDF 预览不可用", "Figure not found",
                                                          "Figure unavailable")) {
  if (!file.exists(html_path)) stop("Narrative report not found: ", html_path)
  html <- paste(readLines(html_path, warn = FALSE), collapse = "\n")

  problems <- character(0)
  has_fig <- grepl("data:image/png;base64,", html, fixed = TRUE)
  n_fig <- if (has_fig) length(gregexpr("data:image/png;base64,", html, fixed = TRUE)[[1]]) else 0
  if (n_fig < min_figures) {
    problems <- c(problems, sprintf("Only %d embedded figures (expected >= %d).", n_fig, min_figures))
  }
  for (marker in missing_markers) {
    if (grepl(marker, html, fixed = TRUE)) {
      problems <- c(problems, paste("Missing-figure marker present in HTML:", marker))
    }
  }
  for (section in required_sections) {
    if (!grepl(section, html, fixed = TRUE)) {
      problems <- c(problems, paste("Required section not found:", section))
    }
  }

  # Failed citer lookups render as "≈ NA", "= NA", "(NA", "（NA" etc. Plain
  # "NA" cannot be grepped (it matches RNA), so anchor on punctuation.
  na_hits <- gregexpr("(≈|＝|=|:|：|\\(|（)\\s*NA([^A-Za-z]|$)", html, perl = TRUE)[[1]]
  n_na <- if (identical(na_hits, -1L)) 0 else length(na_hits)
  if (n_na > 0) {
    message("validate_narrative_report: ", n_na,
            " suspicious inline 'NA' value(s) — inspect whether citer lookups failed ",
            "(a gene/term absent from saved tables renders as NA).")
  }

  if (length(problems)) {
    stop("Narrative report failed validation:\n  - ", paste(problems, collapse = "\n  - "))
  }
  message("Narrative report passed validation: ", n_fig, " embedded figures, ",
          length(required_sections), " required sections present, ",
          n_na, " NA warnings.")
  invisible(list(figures = n_fig, na_warnings = n_na, passed = TRUE))
}

# ---------------------------------------------------------------------------
# Project scaffold
# ---------------------------------------------------------------------------

# Copy the generic narrative-report scaffold into a project's reports/ dir.
# Never overwrites an existing file unless overwrite = TRUE. Returns the path.
scaffold_narrative_report <- function(project_root,
                                      file = file.path(project_root, "reports", "narrative_report.Rmd"),
                                      template = NULL,
                                      overwrite = FALSE) {
  if (!dir.exists(project_root)) stop("Project root not found: ", project_root)
  if (is.null(template)) {
    lib_dir <- Sys.getenv("RNASEQ_LIB_DIR", unset = "")
    if (nzchar(lib_dir) && dir.exists(lib_dir)) {
      template <- file.path(dirname(normalizePath(lib_dir)), "reports", "narrative_report_scaffold.Rmd")
    } else {
      repo_root <- tryCatch(rprojroot::find_root(rprojroot::is_git_root), error = function(e) getwd())
      template <- file.path(repo_root, "reports", "narrative_report_scaffold.Rmd")
    }
  }
  scaffold <- template
  if (!file.exists(scaffold)) stop("Narrative report scaffold not found: ", scaffold)
  if (file.exists(file) && !isTRUE(overwrite)) {
    stop("Refusing to overwrite existing narrative report: ", file,
         " (pass overwrite = TRUE to replace it)")
  }
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  file.copy(scaffold, file, overwrite = overwrite)
  message("Narrative report scaffold written to: ", file, "\n",
          "Edit the project-specific prose (see references/NARRATIVE_REPORT.md), then render with:\n",
          "  rmarkdown::render(\"", basename(file), "\", params = list(project_root = <project>, run_id = <run_id>))")
  invisible(file)
}
