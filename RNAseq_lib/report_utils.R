# Report generation helpers for RNAseq-Templates.
# Renders a self-contained HTML report from the CSV/PDF outputs that the
# notebooks already produce, without re-running any analysis.

# Render the shared Quarto/R Markdown analysis report.
#
# @param outdir Project output directory (contains 1-DEG/, 2-GSEA/, 3-Visualization/).
# @param report_file Path of the HTML report to write.
# @param template Path to the .qmd/.Rmd template (default: reports/analysis_report.qmd
#   resolved relative to the repository root discovered via rprojroot, else CWD).
# @param params Named list of parameters passed to the report (title, author, etc.).
# @return Invisibly, the path to the rendered report.
render_analysis_report <- function(outdir = ".",
                                   report_file = file.path(outdir, "RNAseq_report.html"),
                                   template = NULL,
                                   params = list()) {
  if (is.null(template)) {
    repo_root <- tryCatch(rprojroot::find_root(rprojroot::is_git_root), error = function(e) getwd())
    candidate <- file.path(repo_root, "reports", "analysis_report.qmd")
    template <- if (file.exists(candidate)) candidate else file.path("reports", "analysis_report.qmd")
  }
  if (!file.exists(template)) {
    stop("Report template not found: ", template)
  }

  params <- utils::modifyList(
    list(
      title = "RNA-seq Analysis Report",
      author = "",
      outdir = normalizePath(outdir, mustWork = FALSE)
    ),
    params
  )

  if (requireNamespace("quarto", quietly = TRUE) && nzchar(Sys.which("quarto"))) {
    quarto::quarto_render(
      input = template,
      output_file = basename(report_file),
      execute_params = params,
      quiet = TRUE
    )
    # quarto renders next to the template; move to the requested location.
    produced <- file.path(dirname(template), basename(report_file))
    if (normalizePath(produced, mustWork = FALSE) != normalizePath(report_file, mustWork = FALSE)) {
      file.copy(produced, report_file, overwrite = TRUE)
      unlink(produced)
    }
  } else if (requireNamespace("rmarkdown", quietly = TRUE)) {
    rmarkdown::render(
      input = template,
      output_file = basename(report_file),
      output_dir = dirname(report_file),
      params = params,
      envir = new.env(parent = globalenv()),
      quiet = TRUE
    )
  } else {
    stop("HTML report requires either the 'quarto' CLI (with the R 'quarto' package) ",
         "or the 'rmarkdown' package. Install one of them to enable report generation.")
  }

  if (!file.exists(report_file)) {
    stop("Report rendering did not produce: ", report_file)
  }
  message("HTML report written to: ", report_file)
  invisible(report_file)
}

# List PDF figures under a directory, optionally filtered by a regex pattern.
# Returns paths relative to `base` so the report can embed them with <img>.
list_report_figures <- function(dir, pattern = "\\.pdf$", base = dir) {
  if (!dir.exists(dir)) return(character(0))
  files <- list.files(dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
  if (length(files) == 0) return(character(0))
  # Convert to paths relative to base for portable embedding.
  rel <- sub(paste0("^", normalizePath(base, mustWork = FALSE), .Platform$file.sep, "?"), "",
             normalizePath(files, mustWork = FALSE))
  sort(rel)
}

# Read a CSV if it exists, else return NULL (used to conditionally show tables).
read_csv_safe <- function(path, ...) {
  if (!file.exists(path)) return(NULL)
  tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, ...), error = function(e) NULL)
}
