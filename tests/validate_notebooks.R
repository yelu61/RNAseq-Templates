#!/usr/bin/env Rscript

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required for notebook validation.")
}

repo_root <- tryCatch(
  rprojroot::find_root(rprojroot::is_git_root),
  error = function(e) getwd()
)

notebook_files <- c(
  list.files(file.path(repo_root, "notebooks"), pattern = "[.]ipynb$", full.names = TRUE),
  list.files(file.path(repo_root, "examples"), pattern = "[.]ipynb$", full.names = TRUE, recursive = TRUE)
)

errors <- character(0)
for (file in notebook_files) {
  notebook <- tryCatch(
    jsonlite::fromJSON(file, simplifyVector = FALSE),
    error = function(e) e
  )
  if (inherits(notebook, "error")) {
    errors <- c(errors, paste0(file, ": invalid JSON: ", conditionMessage(notebook)))
    next
  }

  for (i in seq_along(notebook$cells)) {
    cell <- notebook$cells[[i]]
    if (!identical(cell$cell_type, "code")) next
    code <- paste(unlist(cell$source), collapse = "")
    # A previous JSON edit collapsed the whole ThemeEnrichment cell into one
    # comment containing literal "\\n" separators. R parses that as a valid
    # comment, so syntax validation alone cannot catch the lost computation.
    if (grepl("theme_outdir <-", code, fixed = TRUE) &&
        grepl("----\\ntheme_outdir <-", code, fixed = TRUE)) {
      errors <- c(errors, paste0(
        file, " cell ", i,
        ": code contains literal \\n separators; the cell would execute as one comment"
      ))
    }
    parse_error <- tryCatch({
      parse(text = code)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(parse_error)) {
      errors <- c(errors, paste0(file, " cell ", i, ": ", parse_error))
    }
  }
}

library_files <- list.files(file.path(repo_root, "RNAseq_lib"), pattern = "[.]R$", full.names = TRUE)
for (file in library_files) {
  parse_error <- tryCatch({
    parse(file = file)
    NULL
  }, error = function(e) conditionMessage(e))
  if (!is.null(parse_error)) errors <- c(errors, paste0(file, ": ", parse_error))
}

if (length(errors) > 0) stop(paste(errors, collapse = "\n"))
cat("Validated", length(notebook_files), "notebooks and", length(library_files), "R library files.\n")
