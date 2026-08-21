# Run-bundle provenance and registry helpers for RNAseq-Templates.

.run_md5_file <- function(path) {
  if (is.null(path) || length(path) != 1 || is.na(path) || !file.exists(path)) {
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
  files <- unique(normalizePath(files[file.exists(files)], mustWork = TRUE))
  files <- sort(files)
  checksums <- if (length(files)) unname(tools::md5sum(files)) else character(0)
  names(checksums) <- if (length(files)) basename(files) else character(0)
  .run_md5_object(list(revision = revision, file_md5 = checksums))
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
                               backend_files = character(0)) {
  allowed_roles <- c("candidate", "canonical", "sensitivity", "repro_check", "superseded")
  allowed_retention <- c("full", "slim", "metadata_only")
  if (!role %in% allowed_roles) stop("Unsupported run role: ", role)
  if (!retention %in% allowed_retention) stop("Unsupported run retention: ", retention)
  if (!is.list(analysis_config) || is.null(names(analysis_config))) {
    stop("`analysis_config` must be a named list.")
  }
  if (!dir.exists(run_dir)) stop("Run directory not found: ", run_dir)

  run_dir <- normalizePath(run_dir, mustWork = TRUE)
  files <- list.files(run_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[file.info(files)$isdir %in% FALSE]
  run_bytes <- sum(file.info(files)$size, na.rm = TRUE)
  backend <- .run_git_state(backend_root)
  manifest <- data.frame(
    run_id = basename(run_dir),
    status = status,
    role = role,
    parent_run_id = ifelse(is.na(parent_run_id), "", parent_run_id),
    retention = retention,
    completed_at = format(as.POSIXct(completed_at), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    input_file = .run_relative_path(input_file, run_dir),
    input_md5 = .run_md5_file(input_file),
    config_file = .run_relative_path(config_path, run_dir),
    config_md5 = .run_md5_file(config_path),
    analysis_signature = .run_md5_object(analysis_config[sort(names(analysis_config))]),
    backend_revision = backend$revision,
    backend_dirty = backend$dirty,
    backend_signature = .run_backend_signature(backend_files, backend$revision),
    run_bytes = run_bytes,
    change_note = change_note,
    stringsAsFactors = FALSE
  )
  utils::write.csv(manifest, file.path(run_dir, "run_manifest.csv"), row.names = FALSE, na = "")
  invisible(manifest)
}

# Rebuild the project-level registry from per-run manifests. Exact duplicates
# share the input, analysis and backend-code signatures. A full canonical run
# is preferred as owner; other family members point to it through `duplicate_of`.
build_run_registry <- function(runs_dir,
                               path = file.path(runs_dir, "RUN_REGISTRY.csv")) {
  if (!dir.exists(runs_dir)) stop("Runs directory not found: ", runs_dir)
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

  valid_key <- nzchar(registry$input_md5) & nzchar(registry$analysis_signature)
  backend_key <- if ("backend_signature" %in% colnames(registry)) {
    value <- as.character(registry$backend_signature)
    value[is.na(value) | !nzchar(value)] <- "legacy_unknown_backend"
    value
  } else {
    rep("legacy_unknown_backend", nrow(registry))
  }
  keys <- paste(registry$input_md5, registry$analysis_signature, backend_key, sep = "::")
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
