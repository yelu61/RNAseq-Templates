# Run-bundle provenance and registry helpers for RNAseq-Templates.

.run_md5_file <- function(path) {
  if (is.null(path) || length(path) != 1 || is.na(path) || !file.exists(path) || dir.exists(path)) {
    return(NA_character_)
  }
  unname(tools::md5sum(normalizePath(path, mustWork = TRUE))[[1]])
}

.run_canonicalize <- function(x) {
  if (inherits(x, "formula") || is.language(x)) return(paste(deparse(x), collapse = " "))
  if (is.environment(x) || is.function(x)) stop("Run signatures cannot contain environments or functions.")
  if (is.data.frame(x)) {
    out <- lapply(x, .run_canonicalize)
    names(out) <- colnames(x)
    return(list(class = "data.frame", row_names = rownames(x), columns = out))
  }
  if (is.list(x)) {
    out <- lapply(x, .run_canonicalize)
    names(out) <- names(x)
    return(out)
  }
  if (is.factor(x)) return(list(values = as.character(x), levels = levels(x), ordered = is.ordered(x)))
  x
}

.run_md5_object <- function(x) {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(.run_canonicalize(x), tmp, version = 2)
  .run_md5_file(tmp)
}

.run_git_state <- function(root) {
  out <- list(revision = NA_character_, dirty = NA)
  if (is.null(root) || length(root) != 1L || is.na(root) || !dir.exists(root) || !nzchar(Sys.which("git"))) {
    return(out)
  }
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(root)
  revision <- suppressWarnings(tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) character(0)
  ))
  status <- suppressWarnings(tryCatch(
    system2("git", c("status", "--porcelain", "--untracked-files=no"), stdout = TRUE, stderr = FALSE),
    error = function(e) character(0)
  ))
  if (length(revision) && is.null(attr(revision, "status"))) out$revision <- revision[[1]]
  if (!is.na(out$revision)) out$dirty <- length(status) > 0L
  out
}

.run_backend_signature <- function(files = character(0), revision = NA_character_) {
  if (!length(files) || anyNA(vapply(files, .run_md5_file, character(1)))) return(NA_character_)
  files <- unique(normalizePath(files, mustWork = TRUE))
  files <- sort(files)
  checksums <- if (length(files)) unname(tools::md5sum(files)) else character(0)
  names(checksums) <- if (length(files)) basename(files) else character(0)
  .run_md5_object(list(revision = revision, file_md5 = checksums))
}

# This inventory includes the primary matrix plus every declared auxiliary input
# and file-backed reference. Unknown/downloaded references must be declared NA or
# accompanied by inputs_complete = FALSE; omitted references cannot be inferred.
.run_input_inventory <- function(input_file, input_files, run_dir) {
  if (!is.null(input_file) && (!is.character(input_file) || length(input_file) != 1L)) {
    stop("input_file must be one path or NA; use input_files for auxiliary files.")
  }
  if (is.null(input_files)) input_files <- character(0)
  if (!is.character(input_files) || (length(input_files) &&
      (is.null(names(input_files)) || anyNA(names(input_files)) ||
       any(!nzchar(names(input_files))) || anyDuplicated(names(input_files)) ||
       "primary" %in% names(input_files)))) {
    stop("input_files must be a named character vector with unique roles; 'primary' is reserved.")
  }
  paths <- c(primary = if (is.null(input_file)) NA_character_ else unname(input_file), input_files)
  checksums <- vapply(paths, .run_md5_file, character(1))
  specified <- !is.na(paths) & nzchar(paths)
  status <- ifelse(!specified, "unknown", ifelse(!file.exists(paths), "missing",
                   ifelse(is.na(checksums), "unreadable", "present")))
  data.frame(
    role = names(paths),
    path = vapply(paths, .run_relative_path, character(1), run_dir = run_dir),
    md5 = unname(checksums), status = unname(status), stringsAsFactors = FALSE,
    row.names = NULL
  )
}

# Record loaded package versions, including GitHub revisions where available.
# This is a runtime fingerprint, not a restorable environment or a checksum of
# external services, system libraries, or reference data.
.run_runtime_inventory <- function() {
  packages <- sort(loadedNamespaces())
  versions <- vapply(packages, function(package) as.character(getNamespaceVersion(package)), character(1))
  revisions <- vapply(packages, function(package) {
    description <- utils::packageDescription(package)
    revision <- description$RemoteSha
    if (is.null(revision)) revision <- description$GithubSHA1
    if (is.null(revision)) "" else revision
  }, character(1))
  data.frame(
    component = c("R", "platform", paste0("package:", packages)),
    version = c(R.version$version.string, R.version$platform, unname(versions)),
    source_revision = c("", "", unname(revisions)), stringsAsFactors = FALSE
  )
}

# Call before creating output directories or writing the configuration snapshot.
# Input data, config/scripts and an already-open run.log may be present. Reserved
# native result directories are rejected even if the current run disables them.
# output_dirs names explicit output paths (for example a custom WGCNA OUTDIR).
assert_fresh_run_dir <- function(run_dir = ".", output_dirs = character(0)) {
  existing_outputs <- output_dirs[file.exists(output_dirs) | dir.exists(output_dirs)]
  if (length(existing_outputs)) {
    stop("Refusing to reuse an existing analysis output directory: ",
         paste(existing_outputs, collapse = ", "), "\nUse a new run_id.")
  }
  if (!dir.exists(run_dir)) {
    if (file.exists(run_dir)) stop("Run path is not a directory: ", run_dir)
    return(invisible(TRUE))
  }
  entries <- list.files(run_dir, all.files = TRUE, no.. = TRUE)
  blocked <- grepl("^(0-Config$|[1-9][0-9]*-)", entries) |
    entries %in% c("run_manifest.csv", "run_inputs.csv", "run_runtime.csv",
                   "Analysis_summary.txt", "sessionInfo.txt", "RNAseq_report.html",
                   "report_coverage_manifest.csv", "report_review_checklist.csv",
                   "report_validation.csv", "report_interpretation", "report_assets")
  if (any(blocked)) {
    stop("Refusing to reuse a run directory containing analysis outputs: ", run_dir,
         "\nUse a new run_id. Existing artifacts: ", paste(entries[blocked], collapse = ", "))
  }
  invisible(TRUE)
}

.run_relative_path <- function(path, run_dir) {
  if (is.null(path) || length(path) != 1 || is.na(path) || !nzchar(path)) return(NA_character_)
  absolute <- normalizePath(path, mustWork = FALSE)
  run_root <- normalizePath(run_dir, mustWork = FALSE)
  project_root <- if (basename(dirname(run_root)) == "runs") {
    dirname(dirname(dirname(run_root)))
  } else {
    run_root
  }
  prefix <- paste0(project_root, .Platform$file.sep)
  if (startsWith(absolute, prefix)) substring(absolute, nchar(prefix) + 1L) else basename(absolute)
}

initialize_run_lifecycle <- function(envir = parent.frame()) {
  defaults <- list(
    RUN_ROLE = "candidate",
    PARENT_RUN_ID = NA_character_,
    RUN_CHANGE_NOTE = "",
    RUN_RETENTION = "full"
  )
  for (name in names(defaults)) {
    if (!exists(name, envir = envir, inherits = FALSE)) {
      assign(name, defaults[[name]], envir = envir)
    }
  }
  if (!get("RUN_ROLE", envir = envir) %in%
      c("candidate", "canonical", "sensitivity", "repro_check", "superseded")) {
    stop("RUN_ROLE must be candidate/canonical/sensitivity/repro_check/superseded.")
  }
  if (!get("RUN_RETENTION", envir = envir) %in% c("full", "slim", "metadata_only")) {
    stop("RUN_RETENTION must be full/slim/metadata_only.")
  }
  invisible(defaults)
}

write_template_run_manifest <- function(run_dir,
                                        config_path,
                                        input_file,
                                        config_objects,
                                        lib_dir,
                                        runner_file,
                                        envir = parent.frame(),
                                        input_files = NULL,
                                        inputs_complete = !is.null(input_files)) {
  lifecycle_fields <- c("RUN_ROLE", "PARENT_RUN_ID", "RUN_CHANGE_NOTE", "RUN_RETENTION")
  presentation_fields <- c("OUTDIR", "GENERATE_HTML_REPORT", "REPORT_TITLE",
                           "SCAFFOLD_REPORT_INTERPRETATION")
  analysis_objects <- setdiff(config_objects, c(lifecycle_fields, presentation_fields))
  analysis_config <- stats::setNames(lapply(analysis_objects, function(name) {
    if (exists(name, envir = envir, inherits = TRUE)) get(name, envir = envir, inherits = TRUE) else NULL
  }), analysis_objects)
  write_run_manifest(
    run_dir = run_dir,
    config_path = config_path,
    input_file = input_file,
    input_files = input_files,
    inputs_complete = inputs_complete,
    analysis_config = analysis_config,
    role = get("RUN_ROLE", envir = envir, inherits = TRUE),
    parent_run_id = get("PARENT_RUN_ID", envir = envir, inherits = TRUE),
    retention = get("RUN_RETENTION", envir = envir, inherits = TRUE),
    change_note = get("RUN_CHANGE_NOTE", envir = envir, inherits = TRUE),
    backend_root = dirname(normalizePath(lib_dir)),
    backend_files = c(
      if (length(runner_file)) normalizePath(runner_file[[1]], mustWork = TRUE) else character(0),
      list.files(lib_dir, pattern = "[.]R$", full.names = TRUE)
    )
  )
}

# Write one immutable run-level provenance record. This intentionally does not
# update a shared registry: parallel runners can safely write their own
# manifests, and the registry is rebuilt deterministically afterwards.
write_run_manifest <- function(run_dir = ".",
                               config_path,
                               input_file,
                               analysis_config,
                               role = "candidate",
                               parent_run_id = NA_character_,
                               retention = "full",
                               change_note = "",
                               status = "completed",
                               completed_at = Sys.time(),
                               backend_root = NULL,
                               backend_files = character(0),
                               input_files = NULL,
                               inputs_complete = !is.null(input_files)) {
  allowed_roles <- c("candidate", "canonical", "sensitivity", "repro_check", "superseded")
  allowed_retention <- c("full", "slim", "metadata_only")
  if (!role %in% allowed_roles) stop("Unsupported run role: ", role)
  if (!retention %in% allowed_retention) stop("Unsupported run retention: ", retention)
  if (!is.list(analysis_config) || is.null(names(analysis_config))) {
    stop("`analysis_config` must be a named list.")
  }
  if (!dir.exists(run_dir)) stop("Run directory not found: ", run_dir)
  if (!is.logical(inputs_complete) || length(inputs_complete) != 1L || is.na(inputs_complete)) {
    stop("inputs_complete must be TRUE or FALSE.")
  }

  run_dir <- normalizePath(run_dir, mustWork = TRUE)
  provenance_files <- file.path(run_dir, c("run_manifest.csv", "run_inputs.csv", "run_runtime.csv"))
  if (any(file.exists(provenance_files))) stop("Refusing to overwrite existing run provenance in: ", run_dir)
  inputs <- .run_input_inventory(input_file, input_files, run_dir)
  inputs_complete <- inputs_complete && all(inputs$status == "present")
  input_signature <- if (inputs_complete) {
    checksums <- stats::setNames(inputs$md5, inputs$role)
    .run_md5_object(checksums[sort(names(checksums))])
  } else NA_character_
  runtime <- .run_runtime_inventory()
  runtime_signature <- .run_md5_object(runtime)
  backend <- .run_git_state(backend_root)
  backend_signature <- .run_backend_signature(backend_files, backend$revision)
  config_md5 <- .run_md5_file(config_path)
  utils::write.csv(inputs, provenance_files[[2]], row.names = FALSE, na = "")
  utils::write.csv(runtime, provenance_files[[3]], row.names = FALSE, na = "")
  files <- list.files(run_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[file.info(files)$isdir %in% FALSE]
  run_bytes <- sum(file.info(files)$size, na.rm = TRUE)
  manifest <- data.frame(
    manifest_schema_version = 2L,
    run_id = basename(run_dir),
    status = status,
    role = role,
    parent_run_id = ifelse(is.na(parent_run_id), "", parent_run_id),
    retention = retention,
    completed_at = format(as.POSIXct(completed_at), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    input_file = .run_relative_path(input_file, run_dir),
    input_md5 = .run_md5_file(input_file),
    input_manifest_file = "run_inputs.csv",
    input_manifest_md5 = .run_md5_file(provenance_files[[2]]),
    inputs_complete = inputs_complete,
    input_signature = input_signature,
    config_file = .run_relative_path(config_path, run_dir),
    config_md5 = config_md5,
    analysis_signature = .run_md5_object(analysis_config[sort(names(analysis_config))]),
    backend_revision = backend$revision,
    backend_dirty = backend$dirty,
    backend_signature = backend_signature,
    runtime_manifest_file = "run_runtime.csv",
    runtime_signature = runtime_signature,
    runtime_scope = "R_platform_loaded_packages",
    duplicate_eligible = inputs_complete && !is.na(config_md5) &&
      !is.na(backend_signature) && identical(status, "completed"),
    run_bytes = run_bytes,
    change_note = change_note,
    stringsAsFactors = FALSE
  )
  utils::write.csv(manifest, file.path(run_dir, "run_manifest.csv"), row.names = FALSE, na = "")
  invisible(manifest)
}

# Rebuild the registry. Only complete schema-v2 provenance can be grouped by
# declared inputs, analysis, backend code and recorded runtime. Matching these
# signatures is not proof of byte-identical outputs; no pruning is automatic.
build_run_registry <- function(runs_dir,
                               path = file.path(runs_dir, "RUN_REGISTRY.csv")) {
  if (!dir.exists(runs_dir)) stop("Runs directory not found: ", runs_dir)
  runs_dir <- normalizePath(runs_dir, mustWork = TRUE)
  manifests <- list.files(runs_dir, pattern = "^run_manifest[.]csv$",
                          recursive = TRUE, full.names = TRUE)
  manifests <- manifests[dirname(manifests) != normalizePath(runs_dir, mustWork = TRUE)]
  if (!length(manifests)) stop("No run_manifest.csv files found under: ", runs_dir)

  rows <- lapply(manifests, function(file) {
    x <- utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(x) != 1L) stop("Run manifest must contain exactly one row: ", file)
    x
  })
  # Manifest schemas evolve. Fill absent columns before binding so a project
  # can retain old run bundles while adopting newer provenance fields.
  all_columns <- unique(unlist(lapply(rows, colnames), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    for (column in setdiff(all_columns, colnames(x))) x[[column]] <- NA
    x[, all_columns, drop = FALSE]
  })
  registry <- do.call(rbind, rows)
  required <- c("run_id", "completed_at", "input_md5", "analysis_signature")
  missing <- setdiff(required, colnames(registry))
  if (length(missing)) stop("Run manifests are missing columns: ", paste(missing, collapse = ", "))
  registry <- registry[order(registry$completed_at, registry$run_id), , drop = FALSE]
  rownames(registry) <- NULL

  signature_columns <- c("input_md5", "input_signature", "config_md5",
                         "analysis_signature", "backend_signature", "runtime_signature")
  for (column in setdiff(c(signature_columns, "manifest_schema_version", "inputs_complete",
                          "status", "role", "retention"), colnames(registry))) {
    registry[[column]] <- NA
  }
  valid_key <- registry$manifest_schema_version %in% 2L &
    registry$inputs_complete %in% TRUE & registry$status %in% "completed"
  for (column in signature_columns) {
    value <- as.character(registry[[column]])
    valid_key <- valid_key & !is.na(value) & grepl("^[0-9a-f]{32}$", value)
  }
  registry$duplicate_eligible <- valid_key
  keys <- paste(registry$input_signature, registry$analysis_signature,
                registry$backend_signature, registry$runtime_signature, sep = "::")
  registry$duplicate_of <- ""
  role_rank <- match(registry$role,
                     c("canonical", "candidate", "sensitivity", "repro_check", "superseded"))
  role_rank[is.na(role_rank)] <- 99L
  retention_rank <- match(registry$retention, c("full", "slim", "metadata_only"))
  retention_rank[is.na(retention_rank)] <- 99L
  for (key in unique(keys[valid_key])) {
    members <- which(valid_key & keys == key)
    if (length(members) < 2L) next
    owner_order <- order(role_rank[members], retention_rank[members],
                         registry$completed_at[members], registry$run_id[members])
    owner <- members[owner_order[[1]]]
    duplicates <- setdiff(members, owner)
    registry$duplicate_of[duplicates] <- registry$run_id[[owner]]
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(registry, path, row.names = FALSE, na = "")
  invisible(registry)
}
